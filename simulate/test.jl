include("driver-headless.jl")


magnetic_fields = [
    [0,0,0]
    ]

simlist = map(
    b->SimulationParameters("0flux",
    lattice=geom.PyroPrimitive(1,1,1),
    A=[0 0 0 0],
    Jpm=-0.04,
    B=b,
    n_samples=10000
    ),
    magnetic_fields);


for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="tmp",
        sim=sim, 
        integral_params=integration_settings["fast"]
        )
end

