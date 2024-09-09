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

NK = 10000
if length(ARGS) >= 4
    NK = parse(Int64, ARGS[4])
end

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

ip = IntegrationParameters(n_K_samples=NK, integration_method="MC-offset", broaden_factor=2 )

ofname="$(name)%method=$(ip.integration_method)%N=$(ip.n_K_samples)%By=$(bx)at$(point).png"

const Q = SVector{3}(high_symmetry_points[point])
const myB = @SVector [0.0,bx, 0.0]
#const my_Jpm = -0.02
const my_Jpm = -0.02





Egrid = collect(range(0.6,1.4,50));


lats = [
    geom.PyroPrimitive(1,1,1),
    geom.PyroPrimitive(2,1,1),
    geom.PyroPrimitive(2,2,2),
#    geom.PyroPrimitive(1,2,1),
#    geom.PyroPrimitive(1,1,2),
#    geom.PyroPrimitive(3,1,1),
#    geom.PyroPrimitive(5,1,1),
    #geom.PyroPrimitive(7,1,1),
    #geom.PyroPrimitive(1,3,1),
    #geom.PyroPrimitive(1,1,3),
    #geom.PyroPrimitive(2,2,1),
    #geom.PyroPrimitive(2,2,2)
]

zero_A(lat::geom.PyroPrimitive) = zeros(Float64, div(length(lat.tetra_sites), 2), 4)
piflux_A(lat::geom.PyroPrimitive) = construct_landau_gauge(lat, [0 0 0 π; 0 π 0 π; 0 0 0 0])


function randgauge(lat::geom.PyroPrimitive)
    return lattice_gradient(lat,rand(Float64, size(lat.tetra_sites))) .+ 2π*rand(4)'
end

simlist = [
SimulationParameters("$(l.L) rand gauge",
    lattice=l,
    A=zero_A(l), #+ randgauge(l),
    Jpm=my_Jpm,
    B=myB,
    n_samples=1000
    )
    for l in lats
    ]

# for this test's purposes, the lambda used does not need to be correct
# (this value is only good for 0flux)
l0 = 0.12708881513042514
print("Compiling models...  ")
csimlist = [ CompiledModel(s, l0) for s in simlist]
println(collect([calc_lambda(csim) for csim in csimlist]))
print("Done!")



p = plot()
for cs in csimlist
    res = spectral_weight(Q, Egrid, cs, ip,show_progress=true);
    # @assert res.N == ip.n_K_samples

    mean = real.(res.Sqω_pm)./res.N
    var = (res.Sqω_pm2./res.N .- mean.*mean)./res.N

    scale = prod(cs.sim.lat.L)

    plot!(Egrid, mean./scale,
        yerr = sqrt.(var)./scale,
        label=cs.sim.name)
end
print("Saving to ")
println(ofname)

title!("these should overlap")
xlabel!("Energy/Jzz")
ylabel!("intensity")
savefig("testoutput/"*ofname)


