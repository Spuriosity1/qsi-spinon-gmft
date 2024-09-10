include("../src/SpinonStructure.jl")
include("../src/BZMath.jl")


using .BZmath
using .SpinonStructure
using StaticArrays
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using ProgressMeter



println("Testing supercell properties...")
lat_prim=geom.PyroPrimitive(1,1,1)


ch_prim = SpinonStructure.CompiledHamiltonian(sim1)
lat = geom.PyroPrimitive(2,1,1)
