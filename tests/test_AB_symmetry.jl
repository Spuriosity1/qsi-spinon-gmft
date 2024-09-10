include("../src/SpinonStructure.jl")
include("../src/BZMath.jl")


using .BZmath
using .SpinonStructure
using StaticArrays
using SparseArrays
using LinearAlgebra


lat = geom.PyroPrimitive(4,5,6)

function ncells(lat::geom.PyroPrimitive)
    return lat.L[1]*lat.L[2]*lat.L[3]
end


sim = SimulationParameters("0flux",
    lattice= lat,
    A=rand(ncells(lat),4),
    Jpm=-0.001,
    B=[0.,0.,0.],
    n_samples=10000
    )
cmp = SpinonStructure.CompiledHamiltonian(sim)

function assert_sym(M)

    N = ncells(lat)
    
    MA = M[1:N,1:N]
    MB = M[N+1:2N,N+1:2N]
    @assert norm(M[1:N,N+1:2N]) < 1e-15
    @assert norm(M[N+1:2N,1:N]) < 1e-15
    @assert norm(eigvals(MB) - eigvals(MA)) < 1e-15
end


for i=1:1000

    print("Testing... $(i)/1000\r")
    Q = @SVector rand(3)
    assert_sym(Matrix(SpinonStructure.calc_hopping(cmp,Q)))
end

