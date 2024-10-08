module SpinonStructure
"""
Author: Alaric Sanders
Date: 2024-08-14

This module encodes the "physics" part of the code, all other files
    are drivers for confenience.
"""

include("PyrochloreGeometry.jl")
import .PyrochloreGeometry as geom
using StaticArrays
using LinearAlgebra
using Roots
using SparseArrays
using ProgressMeter
using Printf

export calc_fluxes, SimulationParameters, IntegrationParameters, CompiledModel
export geom, calc_lambda, spinon_dispersion, spectral_weight, spectral_weight!
export Sqω_set
export corr_at, broadened_peaks, broadened_peaks!, sim_identifier
export construct_landau_gauge
export lattice_gradient, global_gauge

# this file implements the basic spinon operations.


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


abstract type AbstractCompiledHamiltonian end


"""
    CompiledHamiltonian
Members:
    sim::SimulationParameters

    A sparse matrix of pairs [C, x1-x2] such that the upper diagonal part of H[J1, J2] = C * exp(1.0im*K*(x1-x2)
    coeff_I::Array{Int64}
    coeff_J::Array{Int64}
    coeff_C::Array{ComplexF64}
    coeff_delta::Matrix{Float64}
 
    # the selfconsistent spinon mass
    lambda::Float64
"""
struct CompiledHamiltonian <: AbstractCompiledHamiltonian
    sim::SimulationParameters
    # sparse matrix of pairs [C, delta] such that the upper diagonal part of
    # H[J1, J2](K) = C * exp(1.0im*K*delta)
    coeff_I::Array{Int64}
    coeff_J::Array{Int64}
    coeff_C::Array{ComplexF64}
    coeff_delta::Matrix{Float64}
    # a precomputed set of indices 
    # geom.tetra_idx(sim.lat, sim.lat.tetra_sites[jA] + 2*geom.pyro[mu]) 
    # == nn_index_A[jA, mu]
    nn_index_A::Vector{Vector{Int64}}
    # geom.tetra_idx(sim.lat, sim.lat.tetra_sites[jB] - 2*geom.pyro[mu]) 
    # == nn_index_B[jA, mu]
    nn_index_B::Vector{Vector{Int64}}

    
    function CompiledHamiltonian(_sim::SimulationParameters)
        I = Int64[]
        J = Int64[]
        C = ComplexF64[]
        delta = SVector{3, Float64}[]

        compile_XXZ!(_sim, I, J, C, delta)
        compile_Sx!(_sim, I, J, C, delta)

        nn_idx_A = [
            [geom.tetra_idx(_sim.lat, rA + 2*geom.pyro[mu]) for mu=1:4]
            for rA in geom.A_sites(_sim.lat) 
        ]
        
        nn_idx_B = [
            [geom.tetra_idx(_sim.lat, rA - 2*geom.pyro[mu]) for mu=1:4]
            for rA in geom.B_sites(_sim.lat) 
        ]

        new(_sim,
        I, J, C, reduce(hcat, delta),
        nn_idx_A, nn_idx_B)
    end
end




struct CompiledModel <: AbstractCompiledHamiltonian
    sim::SimulationParameters
    coeff_I::Array{Int64}
    coeff_J::Array{Int64}
    coeff_C::Array{ComplexF64}
    coeff_delta::Matrix{Float64}
    # a precomputed set of indices 
    # geom.tetra_idx(sim.lat, sim.lat.tetra_sites[jA] + 2*geom.pyro[mu]) 
    # == nn_index_A[jA, mu]
    nn_index_A::Vector{Vector{Int64}}
    # geom.tetra_idx(sim.lat, sim.lat.tetra_sites[jB] - 2*geom.pyro[mu]) 
    # == nn_index_B[jA, mu]
    nn_index_B::Vector{Vector{Int64}}

    #the selfconsistent spinon mass
    lambda::Float64
    function CompiledModel(_sim::SimulationParameters)
        ch = CompiledHamiltonian(_sim)
        new(
            _sim,
            ch.coeff_I,
            ch.coeff_J,
            ch.coeff_C,
            ch.coeff_delta,
            ch.nn_index_A,
            ch.nn_index_B,
            calc_lambda(ch)
            )
    end


    function CompiledModel(_sim::SimulationParameters, lambda::Float64)
        ch = CompiledHamiltonian(_sim)
        new(
            _sim,
            ch.coeff_I,
            ch.coeff_J,
            ch.coeff_C,
            ch.coeff_delta,
            ch.nn_index_A,
            ch.nn_index_B,
            lambda
            )
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
#=
@inline function calc_nn_hopping(sim::SimulationParameters, K::Vec3_F64)
    #=
    Calculates the nearest-neighbour hopping matrix elements due to applied magnetic field.
    =#
    # Project B onto the four local z axes
    Bs = map(e-> e[:,3]' * sim.B, geom.axis)
    H = zeros(ComplexF64, length(sim.lat.tetra_sites), length(sim.lat.tetra_sites))
    
    # the nearst neighbour (magnetic) hoppings
    @inbounds begin
    for (J0, x0) in enumerate(geom.A_sites(sim.lat))
        for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)

            hh = Bs[mu]/4 * exp(1.0im*sim.A[J0,mu] - 1.0im*K'*(x1-x0) )
            H[J0, J1] += hh
            H[J1, J0] += conj(hh)
        end
    end
    end
    
    return H
end
=#

#=
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
@inline function calc_xxz_hopping_fast(sim::SimulationParameters, K::Vec3_F64)
    #=
    Calculates the next-neighbour hoppings due to the XXZ like $S^+S^-$ terms, with Jpm = 1.
    =#

    H = zeros(ComplexF64, length(sim.lat.tetra_sites), length(sim.lat.tetra_sites))
    A_sites = geom.A_sites(sim.lat)
    
    @inbounds begin
        
    for (J0, x0) in enumerate(A_sites)
        # the A sites
        for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)
                
            for nu = (mu+1):4
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
        for mu = 1:4
            x1 = x0 - 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)
                
            for nu = mu+1:4
                x2 = x0 - 2*geom.pyro[nu]
                J2 = geom.tetra_idx(sim.lat, x2)
                z= 1/4*exp(1.0im*(-sim.A[J1,mu] + sim.A[J2,nu]) - 1.0im*K'*(x1-x2))
                H[J1, J2] += z
                H[J2, J1] += conj(z)
            end
        end
        
    end
    end

    return H*sim.Jpm
end
=#
#=
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
=#


raw"""
	diagonalise_M(csim, K)
	-> epsilon, U

diagonalises the hopping matrix M, given by the sum of the 
nearest-neighbour and second-neighbour hoppings.
    - K is a k-space point
    - csim is a set of compiled simulation parameters

    Returns:
    - eps, a [M, num_tetra] list of sorted energies
    - U,   a [M, num_tetra, num_tetra] list of eigenvectors such that H[J] @ U[J] =  U[J] @ diag(eps[J])
"""
function diagonalise_M(sim::AbstractCompiledHamiltonian, K::Vec3_F64)
    H = calc_hopping(sim, K) 
    @assert norm(H - H') < 1e-10 
    return eigen!(Matrix(H))
end

    



raw"""
	calc_lambda(sim::AbstractCompiledHamiltonian, nsample, kappa=1)
	-> lambda

Calculates a "chemical potential" lambda that satisfies the averaged self-
consistency condition $\langle \phi^*\phi \rangle = \kappa$.

`bandfunc` is expexted to accept a SVector{Float64,3} and return something 
looking like the output of `eigen!`, i.e. a list of eigenvalues followed by
and eigenvector matrix. (The eigenvectors are not used)

The integral is done using Monte Carlo, samples "nsample" times.
"""
function calc_lambda(csim::AbstractCompiledHamiltonian)
    K = [ geom.primitive_recip_basis * ( 2 .*(@SVector rand(Float64,  3 )) .- 1) for _ in 1:csim.sim.n_samples ]
    
    eps = reduce(vcat, map(k->diagonalise_M(csim, k).values, K)')

    min_lam = -minimum(eps)
    #bandwidth = maximum(eps) + min_lam

    max_lam = csim.sim.kappa^(-2) /2 - minimum(eps)
    
    function constr(λ)
        return sum( (eps .+ λ).^-0.5)/length(eps)/sqrt(2) - csim.sim.kappa
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
	spinon_dispersion(k, sim::SimulationParameters, λ)
	-> E, U

	spinon_dispersion(k, csim::AbstractCompiledHamiltonian)
	-> E, U

Calculates the spinon bands at reciprocal-space point `k` for the given 
parameters. E and U are respectively the spinon energies and the eigenvectors
of the hopping matrix.

k -> either SVector{3,Float64} or Vector{Float64}

Returns U such that 
"""
function spinon_dispersion(k::Vec3_F64, sim::AbstractCompiledHamiltonian, λ::Float64)
    ϵ, U = diagonalise_M(sim, k)
    return sqrt.(2*(ϵ.+ λ )), U
end

function spinon_dispersion(k::Vector{Float64}, sim::AbstractCompiledHamiltonian, λ::Float64)
    ϵ, U = diagonalise_M(sim, SVector{3,Float64}(k))
    return sqrt.(2*(ϵ.+ λ )), U
end

function spinon_dispersion(k::Vec3_F64, sim::CompiledModel)
    ϵ, U = diagonalise_M(sim, k)
    return sqrt.(2*(ϵ.+ sim.lambda )), U
end


function spinon_dispersion(k::Vector{Float64}, sim::CompiledModel)
    ϵ, U = diagonalise_M(sim, SVector{3,Float64}(k))
    return sqrt.(2*(ϵ.+ sim.lambda )), U
end

"""
    IntegrationParameters(n_K_samples::Int)

n_K_samples - number of points to use for the MC integration
"""
@kwdef struct IntegrationParameters
    n_K_samples::Int
    integration_method="MC"
    broaden_factor=1
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
GreenFunction!(G, E, U, i, j)


"""
@inline function GreenFunction!(G, E, U, i, j)
    for l in eachindex(E)
        G[l] = conj(U[i, l]) * U[j,l] / (2*E[l])
    end
end


"""
```
	corr_at(q::Vec3_F64, p::Vec3_F64, sim::AbstractCompiledHamiltonian)
	-> E, Spm, Spp, Smagnetic
```

Calculates the contribution of kspace points(q ± p) to the p integral in 
`corr_Spm`, meaning <S+ S->. Let there be N tetrahedra in `sim`, i.e. N bands.
Returns: 
`E`, an (N, N) matrix of e1 + e2 energies corresponding to Dirac delta peaks
`Spm`, an (N, N) matrix of spectral weights giving the heights of these peaks in <S+(k,w) S-(-k,0)>
`Spp`, an (N, N) matrix of spectral weights giving the heights of these peaks in <S+(k,w) S+(-k,0)>

Notes: 


 I   x1 on the unprimed coords must be conjugated
     relative to x2 on unprimed coords, or else break gauge 
     invariance 
     A \to A + dΓ
     U_{rl} \to e^{iΓ_r} U_{rl}
 II  this conjugation of x1 x2* must be consistent relative to 
     this sign of A
 III The form e^iA_{rA,rAp} e^-iA_{rA+μ,rAp+ν} is certainly correct

"""
function corr_at(Q::Vec3_F64, p1::Vec3_F64, csim::CompiledModel,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing
    )
  
    E1, U1 = spinon_dispersion( p1, csim)
    #p2 = geom.wrap_BZ(csim.sim.lat, p1+Q) # p1+Q
	p2 = p1+Q

    E2, U2 = spinon_dispersion( p2, csim)
    # Both p's must appear with the same sign, or else we break p-> p + delta 
    # invariance (required by gauge symmetry)

	
    QQ_tensor = SMatrix{3,3,Float64}(diagm([1.,1.,1.]) - Q*Q'/(Q'*Q))

	# preallocate everything
    f=length(csim.sim.lat.tetra_sites)
    S_pm = zeros(ComplexF64, f,f)
    S_pp = zeros(ComplexF64, f,f)
    S_magnetic = zeros(Float64, f,f)

    delta_S_pm = zeros(ComplexF64, f,f)
    delta_S_pp = zeros(ComplexF64, f,f)
    delta_S_mag = [ zeros(Float64, f,f) for _ in [1 1; 1 1] ]

    G1 = zeros(ComplexF64, f)
    G2 = zeros(ComplexF64, f)

    G3 = zeros(ComplexF64, f)
    G4 = zeros(ComplexF64, f)

    jB =0
    jpB=0
    
    

    A_sites = geom.A_sites(csim.sim.lat)
    for μ=1:4, ν=1:4
        if g_tensor !== nothing
            R1 = g_tensor * geom.axis[μ]
            R2 = g_tensor * geom.axis[ν]
        end

        for (jA, rA) in enumerate(A_sites), (jpA, rpA) in enumerate(A_sites) 
            jB = csim.nn_index_A[jA][μ]
            jpB = csim.nn_index_A[jpA][ν]
            #@assert jB == geom.tetra_idx(csim.sim.lat, rA + 2*geom.pyro[μ])::Int
            #@assert jpB == geom.tetra_idx(csim.sim.lat, rpA + 2*geom.pyro[ν])::Int

			# <S+ S->
            GreenFunction!(G1, E1, U1, jA, jpA)
			GreenFunction!(G2, E2, U2, jpB, jB)

			#@assert G1 == conj.(U1[jA, :]) .* U1[jpA, :] ./ (2*E1) 
			#@assert G2 == U2[jB, :] .* conj.(U2[jpB, :]) ./ (2*E2)

            for l =1:f, lp=1:f
                delta_S_pm[l,lp] = G1[l]*G2[lp]
            end
            #@assert delta_S_pm == G1*transpose(G2)
 
            delta_S_pm .*= (
				exp(1im*(csim.sim.A[jA,μ]-csim.sim.A[jpA, ν]))
				*exp(1im* ( Q -2*p2)'* (geom.pyro[μ]-geom.pyro[ν])) 
				)
				
			S_pm .+= delta_S_pm

            #################################
			# <S+ S+> 
            GreenFunction!(G3, E1, U1, jA, jpB)
			GreenFunction!(G4, E2, U2, jpA, jB)

            #@assert G3 == conj.(U1[jA, :]) .* U1[jpB, :] ./ (2*E1) 
            #@assert G4 == conj.(U2[jpA, :]) .* U2[jB, :] ./ (2*E2)
 
            for l =1:f, lp=1:f
                delta_S_pp[l,lp] = G3[l]*G4[lp]
            end
            #@assert delta_S_pp == G3*transpose(G4)

            delta_S_pp .*= (
                exp(1im*(Q - 2*p2)'*  geom.pyro[μ] )*
                exp(1im*(-Q - 2*p1)'* geom.pyro[ν]))* 
                #exp(1im*(Q)'* ( geom.pyro[μ] - geom.pyro[ν]))* 
                #exp(-2im*p2'* geom.pyro[μ])*
                #exp(-2im*p1'*geom.pyro[ν]))*
                exp(1im*(csim.sim.A[jA,μ]+csim.sim.A[jpA, ν]))
			S_pp .+= delta_S_pp
            

			if g_tensor !== nothing	
                delta_S_mag[1,1] .=  0.5*real.(delta_S_pp .+ delta_S_pm)
				delta_S_mag[1,2] .=  0.5*imag.(delta_S_pp .- delta_S_pm)
				delta_S_mag[2,1] .=  0.5*imag.(delta_S_pp .+ delta_S_pm)
                delta_S_mag[2,2] .= -0.5*real.(delta_S_pp .- delta_S_pm)

                S_magnetic .+= sum(
                    (R1[:, a]' * QQ_tensor * R2[:, b]) .* delta_S_mag[a, b]
                    for a = 1:2, b = 1:2
                )

			end
			
        end
    end
    E = [e1 + e2 for e1 in E1, e2 in E2]::Matrix{Float64}
    return E, S_pm, S_pp, S_magnetic
end





# puts Lorentzians of weights Snm at energies Enm
function broadened_peaks!(
	Sqω::Vector{Float64},
	Sqω2::Vector{Float64},
	Snm::Matrix{Float64},
	Enm::Matrix{Float64},
	Egrid::Vector{Float64},
	dE::Float64
	)

    for (i,e) in enumerate(Egrid)
        tmp=0.0
        for (E,S) in zip(Enm,Snm)
			tmp += (S*gaussian(e-E, dE))::ComplexF64
		end
        Sqω[i] += tmp
        Sqω2[i] += tmp^2

	end
end


function broadened_peaks!(
	Sqω::Vector{ComplexF64}, 
	Sqω2::Vector{Float64},
	Snm::Matrix{ComplexF64},
	Enm::Matrix{Float64},
	Egrid::Vector{Float64},
	dE::Float64
	)

    for (i,e) in enumerate(Egrid)
        tmp=0.0im
        for (E,S) in zip(Enm,Snm)
			tmp += (S*gaussian(e-E, dE))::ComplexF64
		end
        Sqω[i] += tmp
        Sqω2[i] += real(tmp)^2

	end
end



function broadened_peaks!(
	Sqω::Vector{Float64},
	Snm::Matrix{Float64},
	Enm::Matrix{Float64},
	Egrid::Vector{Float64},
	dE::Float64
	)

	for (E,S) in zip(Enm,Snm)
		for (i,e) in enumerate(Egrid)
            Sqω[i] += S*gaussian(e-E, dE)::Float64
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
            Sqω[i] += S*gaussian(e-E, dE)
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

    # squared vars for extimating variance
    Sqω_pm2::Vector{Float64}
    Sqω_pp2::Vector{Float64}
    Sqω_magnetic2::Vector{Float64}

    Egrid::Vector{Float64}
    bounds::Vector{Float64}
    N::Int

    function Sqω_set(Egrid::Vector{Float64})
        Sqω_pm       = similar(Egrid, ComplexF64)
        Sqω_pp       = similar(Egrid, ComplexF64)
        Sqω_magnetic = similar(Egrid, Float64)        

        Sqω_pm2       = similar(Egrid, Float64)
        Sqω_pp2      = similar(Egrid, Float64)
        Sqω_magnetic2 = similar(Egrid, Float64)        
        new(Sqω_pm, Sqω_pp, Sqω_magnetic,Sqω_pm2,Sqω_pp2,Sqω_magnetic2,
            Egrid,[Inf, -Inf],0)
    end
end

function initzeros!(intensity::Sqω_set)
    intensity.Sqω_pm .= 0
    intensity.Sqω_pp .= 0
    intensity.Sqω_magnetic .= 0


    intensity.Sqω_pm2 .= 0
    intensity.Sqω_pp2 .= 0
    intensity.Sqω_magnetic2 .= 0

    intensity.bounds = [Inf, -Inf]
end


"""
spectral_weight(q, Egrid, sim::CompiledModel,
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
	csim::CompiledModel,
    integral_params::IntegrationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing;
    show_progress=false)
    
    intensity = Sqω_set(Egrid)
 
    spectral_weight!(intensity, 
        q, csim, integral_params, g_tensor,show_progress=show_progress)


    return intensity
end


function get_integration_method(lat::geom.PyroPrimitive, integral_params::IntegrationParameters)
    nsample = nothing
    getp = nothing

    B = geom.reciprocal_basis(lat)

    if integral_params.integration_method == "MC"
        nsample = integral_params.n_K_samples
        getp = (_)-> B * (SVector{3}(rand(3)).-0.5)
	
	elseif integral_params.integration_method == "MC-offset"
        nsample = integral_params.n_K_samples
        getp = (_)-> B * SVector{3}(rand(3))
	elseif integral_params.integration_method == "MC-big"
        nsample = integral_params.n_K_samples
        getp = (_)-> geom.primitive_recip_basis * (SVector{3}(rand(3)).-0.5)
	
    elseif integral_params.integration_method == "grid"
        # round to next largest cube
        nk = Int(ceil(integral_params.n_K_samples^(1/3) ))
        nsample = (nk^3)::Int
        getp = (idx)-> B* ( SVector{3}(
            [ div(idx, nk^2), div(idx, nk) % nk, idx % nk ]
            ) ./nk ) 
    elseif integral_params.integration_method == "grid-big"
        # round to next largest cube
        nk = Int(ceil(integral_params.n_K_samples^(1/3) ))
        nsample = (nk^3)::Int
        getp = (idx)-> geom.primitive_recip_basis* ( SVector{3}(
            [ div(idx, nk^2), div(idx, nk) % nk, idx % nk ]
            ) ./nk ) 
    else
        throw("unrecognised integration method")
    end

    return nsample, getp
end


		# notation: 
		# _pm -> ^{+-}
		# _pp -> ^{++}
		# _rs denotes spinon-site indices
        #
"""
spectral_weight!(
    intensity::Sqω_set, 
    q::Vec3_F64, 
	csim::CompiledModel,
    integral_params::IntegrationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing)
"""
function spectral_weight!(
    intensity::Sqω_set,
    q::Vec3_F64, 
	csim::CompiledModel,
    integral_params::IntegrationParameters,
	g_tensor::Union{Nothing, SMatrix{3,3,Float64}}=nothing;
    show_progress=false)

    initzeros!(intensity) 

    deltaE = maximum(intensity.Egrid[2:end] - intensity.Egrid[1:end-1])

    nsample, getp = get_integration_method(csim.sim.lat, integral_params)
  
    prog = nothing
    if show_progress
        prog=Progress(nsample)
    end

    for idx = 0:(nsample-1)
        p1 =  getp(idx)       # experimental

        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p1, csim, g_tensor)

            broadened_peaks!(intensity.Sqω_pm, intensity.Sqω_pm2,
                    S_pm_rs, E_rs, intensity.Egrid, deltaE*integral_params.broaden_factor )

            intensity.N += 1
            broadened_peaks!(intensity.Sqω_pp, 
                    S_pp_rs, E_rs, intensity.Egrid, deltaE )


            if g_tensor !== nothing
                broadened_peaks!(intensity.Sqω_magnetic, 
                    S_magnetic_rs, E_rs, intensity.Egrid, deltaE )

            end	

            intensity.bounds[1] = min(intensity.bounds[1], reduce(min,  E_rs) ) 
            intensity.bounds[2] = max(intensity.bounds[2], reduce(max,  E_rs) )
        catch e
            if e isa DomainError
                println("Negative dispersion at q=$(q), p1=$(p1)")
                continue
            else
                throw(e)
            end
            
        end


        if show_progress next!(prog) end
    end

    if show_progress finish!(prog) end
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


################################################################################
################################################################################
#  PHYSICS (tight binding) IMPLEMENTATION
#
################################################################################


function compile_XXZ!(sim::SimulationParameters,
    coeff_I::Array{Int64},
    coeff_J::Array{Int64},
    coeff_C::Array{ComplexF64},
    coeff_delta::Array{SVector{3,Float64}}
    )
    
    A_sites = geom.A_sites(sim.lat)
        
    for (J0, x0) in enumerate(A_sites)
        # the A sites
        for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)
                
            for nu = (mu+1):4
                x2 = x0 + 2*geom.pyro[nu]
                J2 = geom.tetra_idx(sim.lat, x2)
                C = sim.Jpm*1/4*exp(1.0im*(-sim.A[J0,mu] + sim.A[J0,nu]))

                push!(coeff_I, J1)
                push!(coeff_J, J2)
                push!(coeff_C, C)
                push!(coeff_delta, -x1+x2)
                # z = 1/4*exp(1.0im*(sim.A[J0,mu] - sim.A[J0,nu]) - 1.0im*K'*(x1-x2))
                @assert x1-x2 == 2*geom.pyro[mu] - 2*geom.pyro[nu]
            end
        end

        # the B sites
        x0 += @SVector [2,2,2]
        J0 += length(A_sites)
        # @assert J0 == tetra_idx(lat, x0) 
        for mu = 1:4
            x1 = x0 - 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)
                
            for nu = mu+1:4
                x2 = x0 - 2*geom.pyro[nu]
                J2 = geom.tetra_idx(sim.lat, x2)
                #z= 1/4*exp(1.0im*(-sim.A[J1,mu] + sim.A[J2,nu]) - 1.0im*K'*(x1-x2))
                #
                C = sim.Jpm*1/4*exp(1.0im*( sim.A[J1,mu] - sim.A[J2,nu]))
                push!(coeff_I, J1)
                push!(coeff_J, J2)

                @assert x1-x2 == -2*geom.pyro[mu] + 2*geom.pyro[nu]
                push!(coeff_C, C) # TEMP TEST
                push!(coeff_delta, -x1+x2)
                #push!(coeff_C, 0 )
            end
        end
		
        
    end

end


function compile_Sx!(sim::SimulationParameters,
    coeff_I::Array{Int64},
    coeff_J::Array{Int64},
    coeff_C::Array{ComplexF64},
    coeff_delta::Array{SVector{3,Float64}}
    )
 
    #=
    Calculates the nearest-neighbour hopping matrix elements due to applied magnetic field.
    =#
    # Project B onto the four local z axes
    Bs = map(e-> e[:,3]' * sim.B, geom.axis)
    
    # the nearst neighbour (magnetic) hoppings
    for (J0, x0) in enumerate(geom.A_sites(sim.lat))
        for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(sim.lat, x1)

            C = Bs[mu]/4 * exp(1.0im*sim.A[J0,mu]) 
            push!(coeff_I, J0)
            push!(coeff_J, J1)
            push!(coeff_C, C)
            push!(coeff_delta, x1-x0)


            #hh = Bs[mu]/4 * exp(1.0im*sim.A[J0,mu] - 1.0im*K'*(x1-x0) )
            #H[J0, J1] += hh
        end
    end
end



@inline function calc_hopping(model::AbstractCompiledHamiltonian, K::Vec3_F64)
    #=
    Calculates the next-neighbour hoppings due to the XXZ like $S^+S^-$ terms, with Jpm = 1.
    =#
    vals_upper = model.coeff_C .* exp.(-1im*model.coeff_delta'*K)::Vector{ComplexF64}
    H = sparse(model.coeff_I, model.coeff_J, vals_upper)
    H += sparse(model.coeff_J, model.coeff_I, conj.(vals_upper))
    return H
end


#############################
#  Convenience functions for gauge transforms
#
"""
global_g must be a row vector
"""
function global_gauge(A_original, global_g)
    return A_original .+ global_g
end

"""
lattice_gradient(lat::geom.PyroPrimitive, G::Vector{Float64})
Calculates the lattice gradient of G, [dG]_{rA,μ}
where rA is an A site position, μ is a sublattice
"""
function lattice_gradient(lat::geom.PyroPrimitive, G::Vector{Float64})
    A_sites = geom.A_sites(lat)
    @assert size(G) == size(lat.tetra_sites)
    retval = zeros(length(A_sites),4)
    for (i,rA) in enumerate(A_sites)
        for μ=1:4
            retval[i,μ] = G[i] - G[geom.tetra_idx(lat, rA + 2 .*geom.pyro[μ])]
        end
    end
    return retval
end

end # end module
