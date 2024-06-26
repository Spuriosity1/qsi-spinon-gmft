include("driver-headless.jl")


Jpm = -0.04

Bmin =sqrt(-9*Jpm/5)
println("minimum B = ",Bmin)

magnetic_fields = [
bb*[0,1,1]/√3 for bb in [Bmin, Bmin + 0.1, Bmin + 0.2]
    ]

lat = geom.PyroPrimitive(2,1,1)

simlist = map(
    b->SimulationParameters("pipi00",
    lattice=lat,
    A=construct_landau_gauge(lat, [0 0 0 π; 0 0 0 0; 0 0 0 0]),
    Jpm=Jpm,
    B=b,
    n_samples = 100000
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

