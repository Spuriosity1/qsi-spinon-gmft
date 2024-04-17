module SpinonStructure



include("PyrochloreGeometry.jl")
import .PyrochloreGeometry as geom
import CSV
using StaticArrays
using LinearAlgebra
using Roots

export load_A, calc_fluxes, SimulationParameters, spinon_dispersion, IntegrationParameters, geom

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
    for (J0, x0) in enumerate(lat.A_sites)
        for mu = 1:4
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
        for mu = 1:4
            x1 = x0 + 2*geom.pyro[mu]
            J1 = geom.tetra_idx(lat, x1)
                
            for nu = (mu+1):4
                x2 = x0 + 2*geom.pyro[nu]
                J2 = geom.tetra_idx(lat, x2)

                z = 1/4*exp(1.0im*(A[J0,mu] - A[J0,nu]) - 1.0im*K'*(x1-x2))
                H[J1, J2] += z
                H[J2, J1] += conj(z)
            end
        end

        # the B sites
        x0 += [2,2,2]
        J0 += length(lat.A_sites)
        # @assert J0 == tetra_idx(lat, x0) 
        for mu = 1:4
            x1 = x0 - 2*geom.pyro[mu]
            J1 = geom.tetra_idx(lat, x1)
                
            for nu = mu+1:4
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
function diagonalise_M(k, A::Matrix{Float64}, Jpm = -1., B = [0.,0.,0.])
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
    K = SVector{3,Float64}(k)

    Jpm = convert(Float64, Jpm)
    B   = convert(Vec3_F64, B)
    
    # infer the lattice size from A
    lat = geom.PyroFCC( round(Int, (size(A)[1]/4)^(1/3)) )
    @assert length(lat.spin_sites) == length(A)
    
    H = Jpm*calc_xxz_hopping_fast(lat, K, A) + calc_nn_hopping(lat, K, A, B)

    @assert norm(H - H') < 1e-10
    
    eps, U = eigen!(Hermitian(H))
    
    
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
    K = rand(Float64, (nsample, 3) )*4π/8
    
    eps = reduce(vcat, map(k->bandf(k)[1],eachrow(K))')

    min_lam = -minimum(eps)
    bandwidth = maximum(eps) + min_lam

    max_lam = kappa^(-2) /2 - minimum(eps)
    
    function constr(λ)
        return mean( (eps .+ λ).^-0.5)/sqrt(2) - kappa
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
      
    function SimulationParameters(x::Dict{String, Any})
        new(x["fluxes"], x["Jpm"], x["B"], x["lambda"], geom.PyroFCC(x["L"]),x["name"])
    end
end

"""
	spinon_dispersion(k, sim)
	-> E, U

Calculates the spinon bands at reciprocal-space point `k` for the given 
parameters. E and U are respectively the spinon energies and the eigenvectors
of the hopping matrix.
"""
function spinon_dispersion(k::Union{Vec3_F64,Vector{Float64}}, sim::SimulationParameters)
    ϵ, U = diagonalise_M(k, sim.A, sim.Jpm, sim.B)    
    return sqrt.(2*(ϵ.+ sim.λ )), U
end



"""
    IntegrationParameters(n_K_samples::Int, BZ_grid_density::Int, broadening_dE::Float64)

n_K_samples - number of points to use for the MC integration
BZ_grid_density - effective length of the system
broadening_DE - lifetime broadening parameter for the Lorentzians
"""
@kwdef struct IntegrationParameters
    n_K_samples::Int
    BZ_grid_density::Int
    broadening_dE::Float64
end

end