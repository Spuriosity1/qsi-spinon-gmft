include("../src/SpinonStructure.jl")
include("../src/BZMath.jl")


using .BZmath
using .SpinonStructure
using StaticArrays
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using Plots

name = ARGS[1]
point = ARGS[2]
bx = parse(Float64, ARGS[3])


const high_symmetry_points = Dict(
    "\\Gamma"=> [0.,0.,0.],
    "X"=> [1.,0.,0.],
    "W"=> [1.,0.5,0.],
    "K"=> [0.75,0.75,0.],
    "L"=> [0.5,0.5,0.5],
    "U"=> [1.0, 0.25,0.25],
    "X2"=> [0.,-1.,0.],
    "W2"=> [0.,1.,0.5],
    "K2"=> [0,0.75,0.75],
    "U2"=> [0.25,1.0, 0.25]
)


const integration_settings = Dict(
    "very_fast" => IntegrationParameters(n_K_samples=10     ),
    "fast" =>   IntegrationParameters(n_K_samples=100       ),
    "slow" =>      IntegrationParameters(n_K_samples=1000   ),
    "very_slow" => IntegrationParameters(n_K_samples=10000  ),
    "ultra_slow" => IntegrationParameters(n_K_samples=100000)
)


ofname="$(name)%By=$(bx)_$(point).png"

const Q = SVector{3}(high_symmetry_points[point])
const myB = @SVector [0,bx, 0]
const my_Jpm = -0.02


ip = integration_settings["very_slow"]


Egrid = collect(range(0,2.2,100));

delta_E = 2.2/100


l221 = geom.PyroPrimitive(2,2,1)
l222 = geom.PyroPrimitive(2,2,2)


A221 =construct_landau_gauge(l221, [0 0 0 π; 0 π 0 π; 0 0 0 0]) + lattice_gradient(l221,rand(Float64, size(l221.tetra_sites))) .+ 2π*rand(4)'
A222 = construct_landau_gauge(l222, [0 0 0 π; 0 π 0 π; 0 0 0 0])+lattice_gradient(l222,rand(Float64, size(l222.tetra_sites))) .+ 2π*rand(4)'

simlist = [
SimulationParameters("221 rand gauge",
    lattice=l221,
    A=A221,
    Jpm=my_Jpm,
    B=myB,
    n_samples=10000
    ),
SimulationParameters("222 loc gauge",
    lattice=l222,
    A=A222,
    Jpm=my_Jpm,
    B=myB,
    n_samples=10000
    )
    ]



print("Compiling models...  ")
csimlist = map(CompiledModel, simlist)
println("Done!")

l0 = csimlist[1].lambda
for j=2:length(csimlist)
    @assert norm(csimlist[j].lambda - l0) < 1e-3
end



p = plot()
for cs in csimlist
    res = spectral_weight(Q, Egrid, cs, ip);

    mean = real.(res.Sqω_pm)./res.N
#    var = res.Sqω_pm2./res.N .- mean.*mean

    scale = prod(cs.sim.lat.L)

    plot!(Egrid, mean./scale,
        #yerr = sqrt.(var)./scale,
        label=cs.sim.name)
end
print("Saving to ")
println(ofname)

title!("these should overlap")
xlabel!("Energy/Jzz")
ylabel!("intensity")
savefig("testoutput/"*ofname)


