include("driver-headless.jl")


Jpm = 0.04

Bmin =sqrt(9*Jpm/5)
println("minimum B = ",Bmin)

magnetic_fields = [
bb*[1,1,0]/√3 for bb in [Bmin, Bmin + 0.1, Bmin + 0.2]
    ]



lat = geom.PyroPrimitive(2,1,1)

simlist = [
    SimulationParameters("pipi00",
    lattice=lat,
    A=construct_landau_gauge(lat, [0 0 0 π; 0 0 0 0; 0 0 0 0]),
    Jpm=Jpm,
    B=b,
    n_samples = 10000
    ) for b in magnetic_fields]

# double check phase is right
for sim in simlist
    g = Jring(Jpm, sim.B )
    phi = calc_fluxes(sim)
    @assert all( g.* cos.(phi) .< 0 )
end

for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["very_slow"]
        )
end

