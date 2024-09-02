include("../src/SimFunctions.jl")

Egrid=collect(range(0,3.2,150))


lat_πflux = geom.PyroPrimitive(2,2,1)
A_piflux = construct_landau_gauge(lat_πflux, [0 0 0 π; 0 π 0 π; 0 0 0 0])

Jpm = 0.05

const Bc = sqrt(9*Jpm/5)

function sim_factory(modB)
    if modB < Bc
        return SimulationParameters("piflux-along",
        A= A_piflux,
        lattice=lat_πflux,
        Jpm=Jpm,
        B=modB*[1.,1.,0.]/√2,
        n_samples=1000,
        kappa=2.0
        )
    else
        return SimulationParameters("0ππ0-along",
        A= [0 π π 0; 0 0 0 0; 0 π π 0; 0 0 0 0 ],
        lattice=geom.PyroPrimitive(1,2,2),
        Jpm=Jpm,
        B=modB*[1.,1.,0.]/√2,
        n_samples=10000,
        kappa=2.0
        )
    
    end
end



G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="output/"

f = integrated_fieldsweep(data_dir, 
                          sim_factory=sim_factory,
                          magnetic_field_strengths=collect(range(0,0.5,50)),
                          #ip=integration_settings["fast"],
                          ip=integration_settings["very_slow"],
                          Egrid=Egrid, g_tensor=G
                         )
println("Saving data to $(f)")





