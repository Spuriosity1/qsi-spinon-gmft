module SpinonStructure



include("PyrochloreGeometry.jl")
import .PyrochloreGeometry as geom
import CSV
using StaticArrays
using LinearAlgebra
using Roots
using Optim
using ProgressMeter

export load_A, calc_fluxes, SimulationParameters, spinon_dispersion, IntegrationParameters, geom, spectral_weight, integrated_specweight, corr_at, broadened_peaks

# The thread switching this does is not desirable at all
LinearAlgebra.BLAS.set_num_threads(1)

# define static convenience types

const Vec3 = geom.Vec3;
const Vec3_F64 = geom.Vec3_F64;

"""
	load_A(filepath::String)

Reads an A file in the plaintext, space-separated format 
x y z A
and returns a ( 4L^3, 4 ) array of gauge connections, where 
A[j, mu] corresponds to the spin at rA_j + geom.pyro[mu]
"""
function load_A(filepath::String)
    file = CSV.File(filepath,ignorerepeated=true, delim=' ');
    
    L = Integer((length(file) / 16)^(1/3))
    
    # Index scheme: A site, sublat
    A = Array{Float64, 2}(undef, 4*L^3, 4)
    
    lat = geom.PyroFCC(L)
    
    for row in file
        psl = row.pyro_sl + 1 # convert to 1-based (PAIN)
        tetra_loc = @SVector[row.X,row.Y,row.Z] .- geom.pyro[psl]
        ti = geom.tetra_idx(lat, tetra_loc)
        A[ti, psl] = row.A
    end
    return A
end

"""
	calc_fluxes(lattice::geom.PyroFCC, A)

Computes the fluxes on `lattice` produces by vector field `A`, in the same 
format as that produced  by `load_A`.
"""
function calc_fluxes(lattice::geom.PyroFCC, A)
    @assert length(A) == length(lattice.spin_sites)
    
    Phi = ones(ComplexF64, length(lattice.A_sites), 4)
    hexa_sites = geom.get_hexagons(lattice)
    
    for (I, row) in enumerate(hexa_sites)
        for (nu, hex) in enumerate(row)
            for (i, (J, sl)) in enumerate(hex)
                Phi[I,nu] *= exp(1im * A[J, sl] * (-1)^i )
            end
        end
    end
   
    return angle.(Phi)
end


raw"""
	calc_nn_hopping(lat, K, A, B)
	->H_B

Calculates the nearest-neigbour part of the hopping matrix $M$ at K-point K, i.e.
the part caused by the (real) magnetic field `B`. Confusingly, `A` refers to _emergent_ magnetic field.
$$
-\sum_{r_A, \mu} \frac{\boldsymbol{\hat{z}}_\mu \cdot \boldsymbol{B}}{4} \left[ \phi_{r_A}^+ e^{iA_{r_A, r_A+b_\mu} } \phi_{r_A + b_\mu}^-  + h.c.\right]
$$
"""
@inline function calc_nn_hopping(lat::geom.PyroFCC, K::Vec3_F64, A::Matrix{Float64}, B::Union{Vector{Float64},Vec3_F64})
    #=
    Calculates the nearest-neighbour hopping matrix elements due to applied magnetic field.
    =#
    # Project B onto the four local z axes
    Bs = map(e-> e[:,3]' * B, geom.axis)
    H = zeros(ComplexF64, length(lat.tetra_sites), length(lat.tetra_sites))
    
    @assert size(A)[1] == length(lat.A_sites)
    # the nearst neighbour (magnetic) hoppings
    @inbounds for (J0, x0) in enumerate(lat.A_sites)
        @inbounds for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(lat, x1)

            hh = Bs[mu]/4 * exp(1.0im*A[J0,mu] - 1.0im*K'*(x1-x0) )
            H[J0, J1] += hh
            H[J1, J0] += conj(hh)
        end
    end
    
    return H
end

raw"""
	calc_xxz_hopping_fast(lat, K, A)
	-> H_xxz

Calculates the second-neighbour hopping caused by the $S^+S^-$ terms at point 
$K$ due to emeergent vector potential $A$.
$$
+ \frac{J_{\pm}}{8}\sum_{r_{\alpha}}\sum_{\mu \neq \nu} 
\phi^{+}_{r_{\alpha}+\mu}\phi^{-}_{r_{\alpha}+\nu}
e^{i\eta_x(A_{r,r+\nu}-A_{r,r+\mu})} 
$$
"""
@inline function calc_xxz_hopping_fast(lat::geom.PyroFCC, K::Vec3_F64, A::Matrix{Float64})
    #=
    Calculates the next-neighbour hoppings due to the XXZ like $S^+S^-$ terms, with Jpm = 1.
    =#

    H = zeros(ComplexF64, length(lat.tetra_sites), length(lat.tetra_sites))
    for (J0, x0) in enumerate(lat.A_sites)
        # the A sites
        @inbounds for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(lat, x1)
                
            @inbounds for nu = (mu+1):4
                x2 = x0 + 2*geom.pyro[nu]
                J2 = geom.tetra_idx(lat, x2)

                z = 1/4*exp(1.0im*(A[J0,mu] - A[J0,nu]) - 1.0im*K'*(x1-x2))
                H[J1, J2] += z
                H[J2, J1] += conj(z)
            end
        end

        # the B sites
        x0 += @SVector [2,2,2]
        J0 += length(lat.A_sites)
        # @assert J0 == tetra_idx(lat, x0) 
        @inbounds for mu = 1:4
            x1 = x0 - 2*geom.pyro[mu]
            J1 = geom.tetra_idx(lat, x1)
                
            @inbounds for nu = mu+1:4
                x2 = x0 - 2*geom.pyro[nu]
                J2 = geom.tetra_idx(lat, x2)
                z= 1/4*exp(1.0im*(-A[J1,mu] + A[J2,nu]) - 1.0im*K'*(x1-x2))
                H[J1, J2] += z
                H[J2, J1] += conj(z)
            end
        end
        
    end

    return H
end


"""
	diagonalise_M(k, A, Jpm, B)
	-> epsilon, U

diagonalises the hopping matrix M, given by the sum of the 
nearest-neighbour and second-neighbour hoppings.
"""
function diagonalise_M(K::Vec3_F64, A::Matrix{Float64}, Jpm = -1., B = [0.,0.,0.])
    #=
    Calculates the spinon hopping matrix and diagonalises it.
    - K is a list of k-space points, e.g. [k1, k2, k3, ... kM]
    - A is a [num_tetra/2, 4] matrix of gauge connections, interpreted as running out of the $A$ sites.
    - Jpm is the XXZ coupling
    - B is the magnetic field in global cubic coords

    Returns:
    - eps, a [M, num_tetra] list of sorted energies
    - U,   a [M, num_tetra, num_tetra] list of eigenvectors such that H[J] @ U[J] =  U[J] @ diag(eps[J])
    =#

    Jpm = convert(Float64, Jpm)
    B   = convert(Vec3_F64, B)
    
    # infer the lattice size from A
    lat = geom.PyroFCC( round(Int, (size(A)[1]/4)^(1/3)) )
    @assert length(lat.spin_sites) == length(A)
    
    H = Jpm*calc_xxz_hopping_fast(lat, K, A) + calc_nn_hopping(lat, K, A, B)

    @assert norm(H - H') < 1e-10
    
    eps, U = eigen!(H)
    
    
    return eps, U
end



raw"""
	calc_lambda(bandfunc, nsample, kappa=1)
	-> lambda

Calculates a "chemical potential" lambda that satisfies the averaged self-
consistency condition $\langle \phi^*\phi \rangle = \kappa$.

`bandfunc` is expexted to accept a SVector{Float64,3} and return something 
looking like the output of `eigen!`, i.e. a list of eigenvalues followed by
and eigenvector matrix. (The eigenvectors are not used)

The integral is done using Monte Carlo, samples "nsample" times.
"""
function calc_lambda(bandf, nsample::Int=1000, kappa::Float64=1)
    K = [ 2π*( (@SVector rand(Float64,  3 )) .-0.5) for _ in 1:nsample ]
    
    eps = reduce(vcat, map(k->bandf(k)[1], K)')

    min_lam = -minimum(eps)
    bandwidth = maximum(eps) + min_lam

    max_lam = kappa^(-2) /2 - minimum(eps)
    
    function constr(λ)
        return sum( (eps .+ λ).^-0.5)/length(eps)/sqrt(2) - kappa
    end

    # function constr_p(λ)
    #     return mean( -0.5 *(eps .+ λ).^-1.5)/sqrt(2)
    # end
    # print(constr(min_lam), constr(max_lam))
    return find_zero(constr, (min_lam, max_lam),Bisection(), naive=true)
    
end

function calc_lambda(A::Matrix{Float64}, Jpm::Float64, 
					 B::Union{Vec3_F64,Vector{Float64}}, 
					 nsample::Int=1000, kappa::Float64=1
					 )

    return calc_lambda(k -> diagonalise_M(k, A, Jpm, B), nsample, kappa)
end


function min_lambda(A::Matrix{Float64}, Jpm::Float64, 
	B::Union{Vec3_F64, Vector{Float64}})
    f = k -> minimum(SpinonStructure.diagonalise_M(k, A, Jpm, B)[1])
	lower = -π/4 .*ones(3)
	upper =  π/4 .*ones(3)
    inner_optimizer = GradientDescent()

    m = Inf
    res = nothing
    for i=1:1000
        r = optimize(f, lower, upper, rand(3)*π/2 .- π/4, Fminbox(inner_optimizer))
        if r.g_converged
            m = min(r.minimum, m)
            res = r
        end
    end
    return res
end

"""
    SimulationParameters
Members:
    A::Matrix{Float64}
    Jpm::Float64
    B::Vec3_F64
    λ::Float64
    lat::geom.PyroFCC

Constructor:
    SimulationParameters(;A, Jpm, B, nsample=10000,  kappa=1.0)
    `A` is the emergent vector potential
    `Jpm` and `B` are Hamiltonian parameters in units of `J_z`
    `lat` is the lattice
    `λ` is the emergent chemical potential.
"""
struct SimulationParameters
    A::Matrix{Float64}
    Jpm::Float64
    B::Vec3_F64
    λ::Float64
    lat::geom.PyroFCC
    name::String

    function SimulationParameters(name; A, Jpm, B, nsample=10000,  kappa=1.0)
        L = round(Int, (size(A)[1]/4)^(1/3))
        new(A, Jpm, B, calc_lambda(A, Jpm, B, nsample, kappa), geom.PyroFCC(L), name )
    end

	function SimulationParameters(name, λ::Float64; A, Jpm, B)
        L = round(Int, (size(A)[1]/4)^(1/3))
		new(A, Jpm, B, λ, geom.PyroFCC(L), name )
	end
      
    function SimulationParameters(x::Dict{String, Any})
        new(x["fluxes"], x["Jpm"], x["B"], x["lambda"], geom.PyroFCC(x["L"]),x["name"])
    end

    function SimulationParameters(other::SimulationParameters, delta_λ::Float64)
        new(other.A, other.Jpm, other.B, other.λ + delta_λ, geom.PyroFCC(other.lat.L), other.name)
    end
end

"""
	spinon_dispersion(k, sim)
	-> E, U

Calculates the spinon bands at reciprocal-space point `k` for the given 
parameters. E and U are respectively the spinon energies and the eigenvectors
of the hopping matrix.
"""
function spinon_dispersion(k::Vec3_F64, sim::SimulationParameters)
    ϵ, U = diagonalise_M(k, sim.A, sim.Jpm, sim.B)
    # if minimum(ϵ.+ sim.λ ) < 0
    #     println("gap closing found near k = $(k)")
    # end
    return sqrt.(2*(ϵ.+ sim.λ )), U
end

function spinon_dispersion(k::Vector{Float64}, sim::SimulationParameters)
    ϵ, U = diagonalise_M(SVector{3,Float64}(k), sim.A, sim.Jpm, sim.B)
    return sqrt.(2*(ϵ.+ sim.λ )), U
end


"""
    IntegrationParameters(n_K_samples::Int, BZ_grid_density::Int, broadening_dE::Float64)

n_K_samples - number of points to use for the MC integration
BZ_grid_density - effective length of the system
broadening_dE - lifetime broadening parameter for the Lorentzians
"""
@kwdef struct IntegrationParameters
    n_K_samples::Int
    broadening_dE::Float64
end


################################################################################
### Calculating the spectral weight
################################################################################

@inline function gaussian(x,σ)
    local N = 1/√(2π)
    return N*exp(-0.5*(x/σ)^2)/σ
end

@inline function Lorentzian(x, Γ)
    local N = 1/π
    Γ /= 2
    return N*1. * Γ/(x^2 + Γ^2)
end


#=
"""
	corr_Spm_at(q::Vec3_F64, p::Vec3_F64, sim::SimulationParameters)
	-> E, S

Calculates the contribution of kspace points(q ± Δ) to the Δ integral in 
`corr_Spm`, meaning <S+ S->. Let there be N tetrahedra in `sim`, i.e. N bands.
Returns: 
`E`, an (N, N) matrix of e1 + e2 energies corresponding to Dirac delta peaks
`S`, an (N, N) matrix of spectral weights giving the heights of these peaks
"""
function _corr_Spm_at(
	E1::Float64, U1::Matrix{ComplexF64},
	E2::Float64, U2::Matrix{ComplexF64},
	q::Vec3_F64, p::Vec3_F64, sim::SimulationParameters)
  
    #E1, U1 = spinon_dispersion( p, sim )
    #E2, U2 = spinon_dispersion( q-p, sim )

    W = zeros(ComplexF64,4,4)
  
    for μ=1:4, ν=1:4
        W[μ,ν] = exp(2im*(q/2 - p)'* (geom.pyro[μ]-geom.pyro[ν]))
    end
  

    local f=length(sim.lat.tetra_sites)
    S = zeros(ComplexF64, f,f)
    for μ=1:4, ν=1:4
        for (jA, rA) in enumerate(sim.lat.A_sites), (jpA, rpA) in enumerate(sim.lat.A_sites)
            jB = geom.tetra_idx(sim.lat, rA + 2*geom.pyro[μ])::Int
            jpB = geom.tetra_idx(sim.lat, rpA + 2*geom.pyro[ν])::Int

            # the "l" bit
            x1 = U1[jA, :] .* conj(U1[jpA, :]) ./ (2*E1) 
            # the "l'" bit
            x2 = conj(U2[jpB, :]) .* U2[jB, :] ./ (2*E2) # <--- calculation says this one
            # x2 = U2[jB, :] .* conj(U2[jpB, :]) ./ (2*E2) # makes no difference 
          
          
            S += W[μ,ν]*exp(1im*(sim.A[jA,μ]-sim.A[jpA, ν])) * x1*x2'
        end
    end
    E = [e1 + e2 for e2 in E2, e1 in E1]::Matrix{Float64}
    return E, S
end




"""
	_corr_Spp_at(
	E1::Float64, U1::Matrix{ComplexF64},
	E2::Float64, U2::Matrix{ComplexF64},
	q::Vec3_F64, p::Vec3_F64, sim::SimulationParameters)
	-> E, S

Calculates the contribution of kspace points(q ± p) to the p integral in 
`corr_Spp`, < S+ S+ >. Let there be N tetrahedra in `sim`, i.e. N bands.
Expects E1, U1 = spinon_dispersion( p, sim )
		E2, U2 = spinon_dispersion( q-p, sim)
Returns: 
`E`, an (N, N) matrix of e1 + e2 energies corresponding to Dirac delta peaks
`S`, an (N, N) matrix of spectral weights giving the heights of these peaks
"""
function _corr_Spp_at(
	E1::Float64, U1::Matrix{ComplexF64},
	E2::Float64, U2::Matrix{ComplexF64},
	q::Vec3_F64, p::Vec3_F64, sim::SimulationParameters)
  
    #E1, U1 = spinon_dispersion( p, sim )
    #E2, U2 = spinon_dispersion( q-p, sim )

    W = zeros(ComplexF64,4,4)
  
    for μ=1:4, ν=1:4
        W[μ,ν] = begin
            exp(2im*(q/2)'* (geom.pyro[μ]-geom.pyro[ν]))*
            exp(-2im*(p)'* (geom.pyro[μ]+geom.pyro[ν]))
        end
    end
  

    local f=length(sim.lat.tetra_sites)
    S = zeros(ComplexF64, f,f)
    for μ=1:4, ν=1:4
        for (jA, rA) in enumerate(sim.lat.A_sites), (jpA, rpA) in enumerate(sim.lat.A_sites)
            jB = geom.tetra_idx(sim.lat, rA + 2*geom.pyro[μ])
            jpB = geom.tetra_idx(sim.lat, rpA + 2*geom.pyro[ν])

            # the "l" bit
            x1 = U1[jA, :] .* conj(U1[jpB, :]) ./ (2*E1) 
            # the "l'" bit
            x2 = U2[jpA, :] .* conj(U2[jB, :]) ./ (2*E2)
          
          
          
            S += W[μ,ν]*exp(1im*(sim.A[jA,μ]+sim.A[jpA, ν])) * x1*x2'
        end
    end
    E = [e1 + e2 for e2 in E2, e1 in E1]::Matrix{Float64}
    return E, S
end
=#



"""
	corr_at(q::Vec3_F64, p::Vec3_F64, sim::SimulationParameters)
	-> E, Spm, Spp, Smagnetic

Calculates the contribution of kspace points(q ± p) to the p integral in 
`corr_Spm`, meaning <S+ S->. Let there be N tetrahedra in `sim`, i.e. N bands.
Returns: 
`E`, an (N, N) matrix of e1 + e2 energies corresponding to Dirac delta peaks
`Spm`, an (N, N) matrix of spectral weights giving the heights of these peaks in <S+(k,w) S-(-k,0)>
`Spp`, an (N, N) matrix of spectral weights giving the heights of these peaks in <S+(k,w) S+(-k,0)>
"""
function corr_at(q::Vec3_F64, p::Vec3_F64, sim::SimulationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing)
  
    E1, U1 = spinon_dispersion( p, sim )
    E2, U2 = spinon_dispersion( q-p, sim )

    local W_pm = zeros(ComplexF64,4,4)
    local W_pp = zeros(ComplexF64,4,4)
  
    @inbounds for μ=1:4, ν=1:4
        W_pm[μ,ν] = exp(2im*(q/2 - p)'* (geom.pyro[μ]-geom.pyro[ν]))

        W_pp[μ,ν] = begin
            exp(2im*(q/2)'* (geom.pyro[μ]-geom.pyro[ν]))*
            exp(-2im*(p)'* (geom.pyro[μ]+geom.pyro[ν]))
        end		
    end

	
	local QQ_tensor = diagm([1.,1.,1.]) - q*q'/(q'*q)	

    local f=length(sim.lat.tetra_sites)
    local S_pm = zeros(ComplexF64, f,f)
    local S_pp = zeros(ComplexF64, f,f)
    local S_magnetic = zeros(Float64, f,f)

    @inbounds for μ=1:4, ν=1:4
        for (jA, rA) in enumerate(sim.lat.A_sites), (jpA, rpA) in enumerate(sim.lat.A_sites)
            jB = geom.tetra_idx(sim.lat, rA + 2*geom.pyro[μ])::Int
            jpB = geom.tetra_idx(sim.lat, rpA + 2*geom.pyro[ν])::Int

			# <S+ S->
            # the "l" bit
            x1 = U1[jA, :] .* conj(U1[jpA, :]) ./ (2*E1) 
            # the "l'" bit
            x2 = conj(U2[jpB, :]) .* U2[jB, :] ./ (2*E2) 
                    
            delta_S_pm = W_pm[μ,ν]*exp(1im*(sim.A[jA,μ]-sim.A[jpA, ν])) * x1*x2'
			S_pm += delta_S_pm
			
			# <S+ S+>
            # the "l" bit
            x3 = U1[jA, :] .* conj(U1[jpB, :]) ./ (2*E1) 
            # the "l'" bit
            x4 = U2[jpA, :] .* conj(U2[jB, :]) ./ (2*E2)

            delta_S_pp = W_pp[μ,ν]*exp(1im*(sim.A[jA,μ]+sim.A[jpA, ν])) * x3*x4'
			S_pp += delta_S_pp


			if g_tensor != nothing

				delta_S_xx =  0.5*real.(delta_S_pp .+ delta_S_pm)
				delta_S_xy =  0.5*imag.(delta_S_pp .- delta_S_pm)
				delta_S_yx =  0.5*imag.(delta_S_pp .+ delta_S_pm)
				delta_S_yy = -0.5*real.(delta_S_pp .- delta_S_pm)

				R1 = g_tensor * geom.axis[μ]
				R2 = g_tensor * geom.axis[ν]

				S_magnetic +=  R1[:,1]' * QQ_tensor * R2[:,1] .* delta_S_xx
				S_magnetic +=  R1[:,1]' * QQ_tensor * R2[:,2] .* delta_S_xy
				S_magnetic +=  R1[:,2]' * QQ_tensor * R2[:,1] .* delta_S_yx
				S_magnetic +=  R1[:,2]' * QQ_tensor * R2[:,2] .* delta_S_yy
	

			end
			
        end
    end
    local E = [e1 + e2 for e2 in E2, e1 in E1]::Matrix{Float64}
    return E, S_pm, S_pp, S_magnetic
end



# puts Lorentzians of weights Snm at energies Enm
function broadened_peaks!(
	Sqω::Union{Vector{ComplexF64}, Vector{Float64}},
	Snm::Union{Matrix{ComplexF64}, Matrix{Float64}},
	Enm::Matrix{Float64},
	Egrid::Vector{Float64},
	dE::Float64
	)

	for (E,S) in zip(Enm,Snm)
		for (i,e) in enumerate(Egrid)
			Sqω[i] += S*Lorentzian(e-E, dE)
		end
	end
end


# puts Lorentzians of weights Snm at energies Enm
# non-mutating version
function broadened_peaks(
	Snm::Union{Matrix{ComplexF64}, Matrix{Float64}},
	Enm::Matrix{Float64},
	Egrid::Vector{Float64},
	dE::Float64
	)
	return reduce( +, [
	[S*Lorentzian(e-E, dE) for e in Egrid]
		for (E,S) in zip(Enm,Snm)
		])		
end



"""
spectral_weight(q, Egrid, sim::SimulationParameters,
integral_params::IntegrationParameters)

Calculates the spectral weight at point `q`.

This performs a Monte Carlo integral of the spin-spin correlators <S+S-> and 
<S+ S+> over the Brillouin zone.

	@param Egrid the energy dgrids to integrate over 
	@param sim the physical parameters of the problem 
	@param integral_params the detials for doing the integration

	@return Sqω_pm The correlator <S^+(q,w) S^-(-q,0)>
	@return Sqω_pp The correlator <S^+(q,w) S^+(-q,0)>
	@return bound, a two-element vector [min,max] bounding the spectral weight

"""
function spectral_weight(q::Vec3_F64, Egrid::Vector{Float64},
	sim::SimulationParameters, integral_params::IntegrationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing)




    # cursed Monte Carlo integration

    bounds = [Inf, -Inf]

	Sqω_pm       = similar(Egrid, ComplexF64)
	Sqω_pp       = similar(Egrid, ComplexF64)
	Sqω_magnetic = similar(Egrid, Float64)

	Sqω_pm       .= 0 
	Sqω_pp       .= 0 
	Sqω_magnetic .= 0 

    
    for _ = 1:integral_params.n_K_samples 
        p = (1 .- 2 .*(@SVector rand(3)))*8π/8         
		#overkill but definitely not too small
		
		# notation: 
		# _pm -> ^{+-}
		# _pp -> ^{++}
		# _rs denotes spinon-site indices
        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, g_tensor)

            broadened_peaks!(Sqω_pm,  S_pm_rs, E_rs, Egrid, integral_params.broadening_dE )
            broadened_peaks!(Sqω_pp,  S_pp_rs, E_rs, Egrid, integral_params.broadening_dE )


            if g_tensor != nothing
                broadened_peaks!(Sqω_magnetic, S_magnetic_rs, E_rs, Egrid, 
                                 integral_params.broadening_dE )

            end	

            bounds[1] = min(bounds[1], reduce(min,  E_rs) ) 
            bounds[2] = max(bounds[2], reduce(max,  E_rs) )
        catch e
            println("Negative dispersion at q=$(q), p=$(p)")
            continue
        end
		# TODO consider doing this in place

    end
    return Sqω_pm, Sqω_pp, Sqω_magnetic, bounds
end


#=
"""
Computes integrated spectral weight over the whole Brillouin zone
"""
function integrated_specweight(sim::SimulationParameters, 
						 integral_params::IntegrationParameters,
        Egrid::Vector{Float64},
        g_tensor::SMatrix{3,3,Float64}
						 )
    Sω_pm = zeros(ComplexF64,size(Egrid))
    Sω_pp = zeros(ComplexF64,size(Egrid))
    Sω_magnetic = zeros(Float64,size(Egrid))

    pr = Progress(integral_params.n_K_samples)
    @Threads.threads for _ = 1:integral_params.n_K_samples  
        q = (1 .- 2 .*(@SVector rand(3)))*4π/8
        p = (1 .- 2 .*(@SVector rand(3)))*4π/8

		E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, g_tensor)

		broadened_peaks!(Sω_pm, S_pm_rs, E_rs, Egrid, integral_params.broadening_dE )
		broadened_peaks!(Sω_pp, S_pp_rs, E_rs, Egrid, integral_params.broadening_dE )

		broadened_peaks!(Sω_magnetic, S_magnetic_rs, E_rs, Egrid,
			integral_params.broadening_dE )


        next!(pr)
    end
    finish!(pr)

    return Sω_pm, Sω_pp, Sω_magnetic
end
=#





end
