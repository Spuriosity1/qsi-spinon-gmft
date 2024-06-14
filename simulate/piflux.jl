include("driver-headless.jl")


magnetic_fields = [
    [0,0,0],
    [1,1,1]*0.1/√3, [0,1,1]*0.1/√2,
    ]

lat = geom.PyroPrimitive(2,2,1)

simlist = map(
    b->SimulationParameters("piflux",
    lattice=lat,
    A=construct_landau_gauge(lat, [0 0 0 π; 0 π 0 π; 0 0 0 0]),
    Jpm=0.3,
    B=b
    ),
    magnetic_fields);


for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["very_slow"]
        )
end

