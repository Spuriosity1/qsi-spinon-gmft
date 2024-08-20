include("../src/SpinonStructure.jl")
include("../src/BZMath.jl")


using .BZmath
using .SpinonStructure
using StaticArrays
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using Plots


const integration_settings = Dict(
    "very_fast" => IntegrationParameters(n_K_samples=10,   broadening_dE=0.1),
    "fast" =>   IntegrationParameters(n_K_samples=100,  broadening_dE=0.05),
    "slow" =>      IntegrationParameters(n_K_samples=1000, broadening_dE=0.02),
    "very_slow" => IntegrationParameters(n_K_samples=10000,broadening_dE=0.02),
    "ultra_slow" => IntegrationParameters(n_K_samples=100000,broadening_dE=0.02)
)


sim0 = SimulationParameters("111 primitive",
    lattice=geom.PyroPrimitive(1,1,1),
    A=zeros(1,4),
    Jpm=-0.046,
    B=[0.,0.,0.],
    n_samples=10000
    )

sim1 = SimulationParameters("221 std_gauge",
    lattice=geom.PyroPrimitive(2,2,1),
    A=zeros(4,4),
    Jpm=-0.046,
    B=[0.,0.,0.],
    n_samples=10000
    )

sim2 = SimulationParameters("221 loc gauge",
    lattice=geom.PyroPrimitive(2,2,1),
    A=[2 2 2 2; -1 -1 -1 -1; 0 0 0 0; 1 1 1 1],
    Jpm=-0.046,
    B=[0.,0.,0.],
    n_samples=10000
    )

sim3 = SimulationParameters("221 NL gauge",
    lattice=geom.PyroPrimitive(2,2,1),
    A=[1 1 0 0; 1 1 0 0; 1 1 0 0; 1 1 0 0],
    Jpm=-0.046,
    B=[0.,0.,0.],
    n_samples=10000
    )

csim0 = CompiledModel(sim0)
csim1 = CompiledModel(sim1)
csim2 = CompiledModel(sim2)
csim3 = CompiledModel(sim3)

@assert norm(csim0.lambda - csim1.lambda) < 1e-3
@assert norm(csim1.lambda - csim2.lambda) < 1e-3
@assert norm(csim1.lambda - csim3.lambda) < 1e-3


Egrid = collect(range(0,2.2,100));
ip = integration_settings["fast"]

# Q = SVector{3}(geom.high_symmetry_points["\\Gamma"])
Q = SVector{3}(geom.high_symmetry_points["X"])



p = plot()
for cs in [csim0, csim1, csim2, csim3]
    res = spectral_weight(Q, Egrid, cs, ip);
    plot!(Egrid, real.(res.Sqω_pm)/prod(cs.sim.lat.L),label=cs.sim.name)
end

name = "comp.png"
if length(ARGS) > 1
    name=ARGS[2]
end
savefig(name)


