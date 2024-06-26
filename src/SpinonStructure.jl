module SpinonStructure

include("PyrochloreGeometry.jl")
import .PyrochloreGeometry as geom
import CSV
using StaticArrays
using LinearAlgebra
using Roots
using Optim
using ProgressMeter
using Printf

export calc_fluxes, SimulationParameters, IntegrationParameters
export geom, calc_lambda, spinon_dispersion, spectral_weight, spectral_weight!
export integrated_specweight, Sqω_set
export corr_at, broadened_peaks, broadened_peaks!, sim_identifier
export construct_landau_gauge


# The thread switching this does is not desirable at all
LinearAlgebra.BLAS.set_num_threads(1)

# define static convenience types

const Vec3 = geom.Vec3;
const Vec3_F64 = geom.Vec3_F64;


"""
    SimulationParameters
Members:
    A::Matrix{Float64}
    Jpm::Float64
    B::Vec3_F64
    lat::geom.PyroGeneric
    name::String

    Contains all basic input parameters needed for a spinon calculation

Constructor:
    SimulationParameters(name;A, Jpm, B, nsample=10000,  kappa=1.0)
    `A` is the emergent vector potential
    `Jpm` and `B` are Hamiltonian parameters in units of `J_z`
    `lat` is the lattice (A is defined with respect to this, see calc_fluxes for details)
"""
struct SimulationParameters
    A::Matrix{Float64}
    Jpm::Float64
    B::Vec3_F64
    lat::geom.PyroGeneric
    kappa::Float64
    n_samples::Int64
    name::String

    function SimulationParameters(name; lattice, A, Jpm, B, kappa=2.0, n_samples=10000)
        @assert length(lattice.spin_sites) == length(A)
        new(A, Jpm, B, lattice, kappa, n_samples, name )
    end

    function SimulationParameters(x::Dict{String, Any})
        lat = begin 
            (x["L"] isa Int) ? 
            geom.PyroFCC(x["L"]) : geom.PyroPrimitive(x["L"]...)
        end
        new(x["gauge"], x["Jpm"], x["B"], lat,
            x["kappa"], x["n_samples"],x["name"])
    end
end


"""
	calc_fluxes(lattice::geom.PyroGeneric, A)

Computes the fluxes on `lattice` produces by vector field `A`. 
A should be a [Lx*Ly*Lz, 4] matrix. A[J,mu] corresponds to 
tetrahedron A_sites(lattice)[J], pyro sublattice mu
"""
function calc_fluxes(lattice::geom.PyroGeneric, A)
    @assert length(A) == length(lattice.spin_sites)
    
    Phi = ones(ComplexF64, length( geom.A_sites(lattice)), 4)
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


function calc_fluxes(sim::SimulationParameters)
    return calc_fluxes(sim.lat, sim.A)
end


raw"""
	calc_nn_hopping(sim::SimulationParameters, K::Vec3_F64)
	->H_B

Calculates the nearest-neigbour part of the hopping matrix $M$ at K-point K, i.e.
the part caused by the (real) magnetic field `B`. Confusingly, `A` refers to _emergent_ magnetic field.
$$
-\sum_{r_A, \mu} \frac{\boldsymbol{\hat{z}}_\mu \cdot \boldsymbol{B}}{4} \left[ \phi_{r_A}^+ e^{iA_{r_A, r_A+b_\mu} } \phi_{r_A + b_\mu}^-  + h.c.\right]
$$
"""
# @inline function calc_nn_hopping(lat::geom.PyroGeneric, K::Vec3_F64, A::Matrix{Float64}, B::Union{Vector{Float64},Vec3_F64})
@inline function calc_nn_hopping(sim::SimulationParameters, K::Vec3_F64)
    #=
    Calculates the nearest-neighbour hopping matrix elements due to applied magnetic field.
    =#
    # Project B onto the four local z axes
    Bs = map(e-> e[:,3]' * sim.B, geom.axis)
    H = zeros(ComplexF64, length(sim.lat.tetra_sites), length(sim.lat.tetra_sites))
    
    # the nearst neighbour (magnetic) hoppings
    @inbounds for (J0, x0) in enumerate(geom.A_sites(sim.lat))
        @inbounds for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)

            hh = Bs[mu]/4 * exp(1.0im*sim.A[J0,mu] - 1.0im*K'*(x1-x0) )
            H[J0, J1] += hh
            H[J1, J0] += conj(hh)
        end
    end
    
    return H
end


raw"""
	calc_xxz_hopping_fast(sim , K)
	-> H_xxz

Calculates the second-neighbour hopping caused by the $S^+S^-$ terms at point 
$K$ due to emeergent vector potential $A$.
$$
+ \frac{J_{\pm}}{8}\sum_{r_{\alpha}}\sum_{\mu \neq \nu} 
\phi^{+}_{r_{\alpha}+\mu}\phi^{-}_{r_{\alpha}+\nu}
e^{i\eta_x(A_{r,r+\nu}-A_{r,r+\mu})} 
$$
"""
# @inline function calc_xxz_hopping_fast(lat::geom.PyroGeneric, K::Vec3_F64, A::Matrix{Float64})
@inline function calc_xxz_hopping_fast(sim::SimulationParameters, K::Vec3_F64)
    #=
    Calculates the next-neighbour hoppings due to the XXZ like $S^+S^-$ terms, with Jpm = 1.
    =#

    H = zeros(ComplexF64, length(sim.lat.tetra_sites), length(sim.lat.tetra_sites))
    A_sites = geom.A_sites(sim.lat)
    
    for (J0, x0) in enumerate(A_sites)
        # the A sites
        @inbounds for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)
                
            @inbounds for nu = (mu+1):4
                x2 = x0 + 2*geom.pyro[nu]
                J2 = geom.tetra_idx(sim.lat, x2)

                z = 1/4*exp(1.0im*(sim.A[J0,mu] - sim.A[J0,nu]) - 1.0im*K'*(x1-x2))
                H[J1, J2] += z
                H[J2, J1] += conj(z)
            end
        end

        # the B sites
        x0 += @SVector [2,2,2]
        J0 += length(A_sites)
        # @assert J0 == tetra_idx(lat, x0) 
        @inbounds for mu = 1:4
            x1 = x0 - 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)
                
            @inbounds for nu = mu+1:4
                x2 = x0 - 2*geom.pyro[nu]
                J2 = geom.tetra_idx(sim.lat, x2)
                z= 1/4*exp(1.0im*(-sim.A[J1,mu] + sim.A[J2,nu]) - 1.0im*K'*(x1-x2))
                H[J1, J2] += z
                H[J2, J1] += conj(z)
            end
        end
        
    end

    return H*sim.Jpm
end


raw"""
	diagonalise_M(sim, K)
	-> epsilon, U

diagonalises the hopping matrix M, given by the sum of the 
nearest-neighbour and second-neighbour hoppings.
    - K is a k-space point
    - sim is a set of simulatipn parameters

    Returns:
    - eps, a [M, num_tetra] list of sorted energies
    - U,   a [M, num_tetra, num_tetra] list of eigenvectors such that H[J] @ U[J] =  U[J] @ diag(eps[J])
"""
function diagonalise_M(sim::SimulationParameters, K::Vec3_F64)

    
    H = calc_xxz_hopping_fast(sim, K) + calc_nn_hopping(sim, K)

    @assert norm(H - H') < 1e-10
    
    eps, U = eigen!(H)
    
    
    return eps, U
end



raw"""
	calc_lambda(sim::SimulationParameters, nsample, kappa=1)
	-> lambda

Calculates a "chemical potential" lambda that satisfies the averaged self-
consistency condition $\langle \phi^*\phi \rangle = \kappa$.

`bandfunc` is expexted to accept a SVector{Float64,3} and return something 
looking like the output of `eigen!`, i.e. a list of eigenvalues followed by
and eigenvector matrix. (The eigenvectors are not used)

The integral is done using Monte Carlo, samples "nsample" times.
"""
function calc_lambda(sim::SimulationParameters)
    K = [ 2π*( (@SVector rand(Float64,  3 )) .-0.5) for _ in 1:sim.n_samples ]
    
    eps = reduce(vcat, map(k->diagonalise_M(sim, k)[1], K)')

    min_lam = -minimum(eps)
    bandwidth = maximum(eps) + min_lam

    max_lam = sim.kappa^(-2) /2 - minimum(eps)
    
    function constr(λ)
        return sum( (eps .+ λ).^-0.5)/length(eps)/sqrt(2) - sim.kappa
    end

    # function constr_p(λ)
    #     return mean( -0.5 *(eps .+ λ).^-1.5)/sqrt(2)
    # end
    # print(constr(min_lam), constr(max_lam))
    return find_zero(constr, (min_lam, max_lam),Bisection(), naive=true)
end


# function min_lambda(A::Matrix{Float64}, Jpm::Float64, 
# 	B::Union{Vec3_F64, Vector{Float64}}, n_tries=1000)
#     f = k -> minimum(SpinonStructure.diagonalise_M(SVector{3}(k), A, Jpm, B)[1])
# 	lower = -π/4 .*ones(3)
# 	upper =  π/4 .*ones(3)
#     inner_optimizer = GradientDescent()

#     m = Inf
#     res = nothing
#     prog = Progress(n_tries)
#     @Threads.threads for i=1:n_tries
#         r = optimize(f, lower, upper, rand(3)*π/2 .- π/4, Fminbox(inner_optimizer))
#         if r.g_converged
#             m = min(r.minimum, m)
#             res = r
#         end
#         next!(prog)
#     end
#     finish!(prog)
#     return res
# end





"""
A 'hash function' providing an informative filename for simulation data
"""
function sim_identifier(sim::SimulationParameters)
    return @sprintf("?name=%s?J_pm=%.3f?B=[%.3f,%.3f,%.3f]",
        sim.name,sim.Jpm,sim.B[1],sim.B[2],sim.B[3]) 
end

"""
	spinon_dispersion(k, sim, λ)
	-> E, U

Calculates the spinon bands at reciprocal-space point `k` for the given 
parameters. E and U are respectively the spinon energies and the eigenvectors
of the hopping matrix.
"""
function spinon_dispersion(k::Vec3_F64, sim::SimulationParameters, λ::Float64)
    ϵ, U = diagonalise_M(sim, k)
    # if minimum(ϵ.+ sim.λ ) < 0
    #     println("gap closing found near k = $(k)")
    # end
    return sqrt.(2*(ϵ.+ λ )), U
end

function spinon_dispersion(k::Vector{Float64}, sim::SimulationParameters, λ::Float64)
    ϵ, U = diagonalise_M(sim, SVector{3,Float64}(k))
    return sqrt.(2*(ϵ.+ λ )), U
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
        λ::Float64,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing)
  
    E1, U1 = spinon_dispersion( p, sim, λ)
    E2, U2 = spinon_dispersion( q-p, sim, λ)

	
    QQ_tensor = SMatrix{3,3,Float64}(diagm([1.,1.,1.]) - q*q'/(q'*q))

    f=length(sim.lat.tetra_sites)
    S_pm = zeros(ComplexF64, f,f)
    S_pp = zeros(ComplexF64, f,f)
    S_magnetic = zeros(Float64, f,f)

    jB =0
    jpB=0
    

    delta_S_pm = Array{ComplexF64}(undef, f,f)
    delta_S_pp = Array{ComplexF64}(undef, f,f)
    

    A_sites = geom.A_sites(sim.lat)
    @inbounds for μ=1:4, ν=1:4
        for (jA, rA) in enumerate(A_sites), (jpA, rpA) in enumerate(A_sites)
            jB = geom.tetra_idx(sim.lat, rA + 2*geom.pyro[μ])::Int
            jpB = geom.tetra_idx(sim.lat, rpA + 2*geom.pyro[ν])::Int

			# <S+ S->
            # the "l" bit
            x1 = U1[jA, :] .* conj(U1[jpA, :]) ./ (2*E1) 
            # the "l'" bit
            x2 = conj(U2[jpB, :]) .* U2[jB, :] ./ (2*E2) 
                    
            delta_S_pm = (
                exp(2im*(q/2 - p)'* (geom.pyro[μ]-geom.pyro[ν]))
                )*exp(1im*(sim.A[jA,μ]-sim.A[jpA, ν])) * x1*x2'
			S_pm .+= delta_S_pm
			
			# <S+ S+>
            # the "l" bit
            x3 = U1[jA, :] .* conj(U1[jpB, :]) ./ (2*E1) 
            # the "l'" bit
            x4 = U2[jpA, :] .* conj(U2[jB, :]) ./ (2*E2)

            delta_S_pp = (
                exp(2im*(q/2)'* (geom.pyro[μ]-geom.pyro[ν]))*
                exp(-2im*(p)'* (geom.pyro[μ]+geom.pyro[ν]))
                )*exp(1im*(sim.A[jA,μ]+sim.A[jpA, ν])) * x3*x4'
			S_pp .+= delta_S_pp


			if g_tensor != nothing

				delta_S_xx =  0.5*real.(delta_S_pp .+ delta_S_pm)
				delta_S_xy =  0.5*imag.(delta_S_pp .- delta_S_pm)
				delta_S_yx =  0.5*imag.(delta_S_pp .+ delta_S_pm)
				delta_S_yy = -0.5*real.(delta_S_pp .- delta_S_pm)

				R1 = g_tensor * geom.axis[μ]
				R2 = g_tensor * geom.axis[ν]

				S_magnetic .+=  R1[:,1]' * QQ_tensor * R2[:,1] .* delta_S_xx
				S_magnetic .+=  R1[:,1]' * QQ_tensor * R2[:,2] .* delta_S_xy
				S_magnetic .+=  R1[:,2]' * QQ_tensor * R2[:,1] .* delta_S_yx
				S_magnetic .+=  R1[:,2]' * QQ_tensor * R2[:,2] .* delta_S_yy
	

			end
			
        end
    end
    E = [e1 + e2 for e2 in E2, e1 in E1]::Matrix{Float64}
    return E, S_pm, S_pp, S_magnetic
end





# puts Lorentzians of weights Snm at energies Enm
function broadened_peaks!(
	Sqω::Vector{Float64},
	Snm::Matrix{Float64},
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


function broadened_peaks!(
	Sqω::Vector{ComplexF64},
	Snm::Matrix{ComplexF64},
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

mutable struct Sqω_set
    Sqω_pm::Vector{ComplexF64}
    Sqω_pp::Vector{ComplexF64}
    Sqω_magnetic::Vector{Float64}
    Egrid::Vector{Float64}
    bounds::Vector{Float64}

    function Sqω_set(Egrid::Vector{Float64})
        Sqω_pm       = similar(Egrid, ComplexF64)
        Sqω_pp       = similar(Egrid, ComplexF64)
        Sqω_magnetic = similar(Egrid, Float64)        

        new(Sqω_pm, Sqω_pp, Sqω_magnetic, Egrid,[Inf, -Inf])
    end
end

function initzeros!(intensity::Sqω_set)
    intensity.Sqω_pm .= 0
    intensity.Sqω_pp .= 0
    intensity.Sqω_magnetic .= 0

    intensity.bounds = [Inf, -Inf]
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
	sim::SimulationParameters, λ::Float64,
    integral_params::IntegrationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing)
    
    intensity = Sqω_set(Egrid)
 
    spectral_weight!(intensity, 
        q, sim, λ, integral_params, g_tensor)


    return intensity
end

function spectral_weight!(
    intensity::Sqω_set,
    q::Vec3_F64, 
	sim::SimulationParameters, λ::Float64,
    integral_params::IntegrationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing)
    # cursed Monte Carlo integration

    initzeros!(intensity) 

    for _ = 1:integral_params.n_K_samples 
        p = (1 .- 2 .*(@SVector rand(3)))*2π
		#overkill but definitely not too small
		
		# notation: 
		# _pm -> ^{+-}
		# _pp -> ^{++}
		# _rs denotes spinon-site indices
        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, λ, g_tensor)

            broadened_peaks!(intensity.Sqω_pm,  S_pm_rs, E_rs, intensity.Egrid, integral_params.broadening_dE )
            broadened_peaks!(intensity.Sqω_pp,  S_pp_rs, E_rs, intensity.Egrid, integral_params.broadening_dE )


            if g_tensor != nothing
                broadened_peaks!(intensity.Sqω_magnetic, S_magnetic_rs, E_rs, intensity.Egrid, 
                                 integral_params.broadening_dE )

            end	

            intensity.bounds[1] = min(intensity.bounds[1], reduce(min,  E_rs) ) 
            intensity.bounds[2] = max(intensity.bounds[2], reduce(max,  E_rs) )
        catch e
            if e isa DomainError
                println("Negative dispersion at q=$(q), p=$(p)")
                continue
            else
                throw(e)
            end
            
        end

    end
end

"""
function construct_landau_gauge(lattice::geom.PyroPrimitive, α)
Input: primitive lattice, and a field gradient α
α should be a 3 x 4 matrix such that 
A_μ (I_i e_i ) = I_i α_iμ
Output: A, a (V,4) matrix contianing the gauge fields
"""
function construct_landau_gauge(lattice::geom.PyroPrimitive, α)
    A_tet_sites = geom.A_sites(lattice)
    
    A = zeros(Float64, length(A_tet_sites),4)
    for (J, tA) in enumerate(A_tet_sites)
        I, _ = geom.tetra_IDX(lattice, tA) # The three-tuple that does the thing
        A[J,:] = I' * α
    end
    # consistency
    # println(exp.(1im* lat.L'*α))
    @assert all(abs.(exp.(1im* lattice.L'*α) .- 1.) .< 1e-10) """
        phases do not match when translating across the unit cell boundary,
        all entries in α[:,i] must be integer multiples of 2π/lat.L[i]
        """
        
    return A
end



end # end module
