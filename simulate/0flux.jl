import Pkg; Pkg.activate(".")
include("driver-headless.jl")


magnetic_fields = [
    [0,0,0],
    [1,1,1]*0.1/√3, [1,1,0]*0.1/√2,
    ]

    Jpm = -0.04

    # double check phase is right
for B in magnetic_fields
    @assert all(Jring(-0.04, B).< 0)
end


simlist = [
    SimulationParameters("0flux",
    lattice=geom.PyroPrimitive(1,1,1),
    A=[0 0 0 0],
    Jpm=-0.04,
    B=b,
    n_samples=10000
    ) for b in magnetic_fields ]


for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["ultra_slow"]
        )
end

