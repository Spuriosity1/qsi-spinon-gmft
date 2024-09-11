include("../src/SpinonStructure.jl")
include("../src/BZMath.jl")


using .BZmath
using .SpinonStructure
using LaTeXStrings
using StaticArrays
using SparseArrays
using LinearAlgebra
using Plots

if length(ARGS) < 3
    println("USAGE: name point bx [NK]")
    exit(1)
end

name = ARGS[1]
point = ARGS[2]
bx = parse(Float64, ARGS[3])

NK = 1000
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

ip = IntegrationParameters(n_K_samples=NK, integration_method="grid", broaden_factor=2 )

ofname="$(name)%method=$(ip.integration_method)%N=$(ip.n_K_samples)%By=$(bx)at$(point)"

const Q = SVector{3}(high_symmetry_points[point])
const myB = @SVector [-bx*0.5,bx, bx*0.3]
const my_Jpm = -0.02





Egrid = collect(range(0.6,1.4,50));

lat = geom.PyroPrimitive(1,2,3)
#lat = geom.PyroPrimitive(1,1,1)
#lat = geom.PyroPrimitive(2,3,2)

zero_A(lat::geom.PyroPrimitive) = zeros(Float64, div(length(lat.tetra_sites), 2), 4)
piflux_A(lat::geom.PyroPrimitive) = construct_landau_gauge(lat, [0 0 0 π; 0 π 0 π; 0 0 0 0])


function randgauge(lat::geom.PyroPrimitive)
    return lattice_gradient(lat,rand(Float64, size(lat.tetra_sites))) 
end

simlist = [
SimulationParameters("reference",
    lattice=lat,
    A=zero_A(lat),
    Jpm=my_Jpm,
    B=myB,
    n_samples=50
    ),
SimulationParameters("random div free",
    lattice=lat,
    A=zero_A(lat)+ randgauge(lat),
    Jpm=my_Jpm,
    B=myB,
    n_samples=50
    ),
SimulationParameters("random shift",
    lattice=lat,
    A=zero_A(lat).+ 2π*rand(4)',
    Jpm=my_Jpm,
    B=myB,
    n_samples=50
    )
    ]

# for this test's purposes, the lambda used does not need to be correct
# (this value is only good for 0flux)
l0 = 0.18
print("Compiling models...  ")
csimlist = [ CompiledModel(s, l0) for s in simlist]
println(collect(cs.lambda for cs in csimlist))
print("Done!")

results = [spectral_weight(Q, Egrid, cs, ip,show_progress=true) 
    for cs in csimlist ]


p = plot()
for (cs,res) in zip(csimlist, results)

    mean = real.(res.Sqω_pm)./res.N
    var = (res.Sqω_pm2./res.N .- mean.*mean)./res.N

    scale = prod(cs.sim.lat.L)

    plot!(Egrid, mean./scale,
        yerr = sqrt.(var)./scale,
        label=cs.sim.name)
end
println("Saving to $(ofname)_pm_re.png")

title!(L"\rm{Re}\langle S^+S^-\rangle")
xlabel!("Energy/Jzz")
ylabel!("intensity")
savefig("testoutput/$(ofname)_pm_re.png")


p = plot()
for (cs,res) in zip(csimlist, results)

    mean = imag.(res.Sqω_pm)./res.N
    var = (res.Sqω_pm2./res.N .- mean.*mean)./res.N

    scale = prod(cs.sim.lat.L)

    plot!(Egrid, mean./scale,
        yerr = sqrt.(var)./scale,
        label=cs.sim.name)
end
println("Saving to $(ofname)_pm_im.png")

title!(L"\rm{Im}\langle S^+S^-\rangle")
xlabel!("Energy/Jzz")
ylabel!("intensity")
savefig("testoutput/$(ofname)_pm_im.png")

p = plot()
for (cs,res) in zip(csimlist, results)

    mean = real.(res.Sqω_pp)./res.N
    var = [max(0, v) for v in (res.Sqω_pp2./res.N .- mean.*mean)./res.N]

    scale = prod(cs.sim.lat.L)

    plot!(Egrid, mean./scale,
        yerr = sqrt.(var)./scale,
        label=cs.sim.name)
end
println("Saving to $(ofname)_pp_Re.png")

title!(L"\rm{Re}\langle S^+S^+\rangle")
xlabel!("Energy/Jzz")
ylabel!("intensity")
savefig("testoutput/$(ofname)_pp_Re.png")



p = plot()
for (cs,res) in zip(csimlist, results)

    mean = imag.(res.Sqω_pp)./res.N
    var = [max(0, v) for v in (res.Sqω_pp2./res.N .- mean.*mean)./res.N]

    scale = prod(cs.sim.lat.L)

    plot!(Egrid, mean./scale,
        yerr = sqrt.(var)./scale,
        label=cs.sim.name)
end
println("Saving to $(ofname)_pp_Im.png")


title!(L"\rm{Im}\langle S^+S^+\rangle")
xlabel!("Energy/Jzz")
ylabel!("intensity")
savefig("testoutput/$(ofname)_pp_Im.png")

