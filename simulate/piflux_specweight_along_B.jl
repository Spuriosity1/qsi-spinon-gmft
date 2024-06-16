include("../src/SimFunctions.jl")

Egrid=collect(range(0,3.2,150))


lat = geom.PyroPrimitive(2,2,1)
A_piflux = construct_landau_gauge(lat, [0 0 0 π; 0 π 0 π; 0 0 0 0])

function sim_factory(modB)
    return SimulationParameters("piflux-along",
    A= A_piflux,
    lattice=lat,
    Jpm=0.3,
    B=modB*[1.,1.,1.]/√3,
    nsample=10000,
    kappa=2.0
    )
end


G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="output/"

f = integrated_fieldsweep(data_dir, 
                          sim_factory=sim_factory,
                          magnetic_field_strengths=collect(range(0,0.5,50)),
                          ip=integration_settings["very_slow"],
                          Egrid=Egrid, g_tensor=G
                         )
println("Saving data to $(f)")





