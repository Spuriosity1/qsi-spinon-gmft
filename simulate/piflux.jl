include("driver-headless.jl")


magnetic_fields = [
    [0,0,0],
    [1,1,1]*0.1/√3, [0,1,1]*0.1/√2,
    ]

const Jpm = 0.3
const Bc = sqrt(9*Jpm/5)

@assert all( norm(b) < Bc for b in magnetic_fields)

lat = geom.PyroPrimitive(2,2,1)

simlist = [
    SimulationParameters("piflux",
    lattice=lat,
    A=construct_landau_gauge(lat, [0 0 0 π; 0 π 0 π; 0 0 0 0]),
    Jpm=Jpm,
    B=b
    ) for b in
    magnetic_fields]


# double check phase is right
for sim in simlist
    g = Jring(Jpm, sim.B )
    phi = calc_fluxes(sim)
    @assert all( g.* cos.(phi) < 0 )
end


for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["slow"],
        k_density_spinon_dispersion=120,
        calc_specweight=true,
        calc_integrated=true
        )
end

