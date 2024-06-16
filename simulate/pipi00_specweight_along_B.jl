include("../src/SimFunctions.jl")

Egrid=collect(range(0,3.2,150))


Jpm = -0.05

Bmin =sqrt(-9*Jpm/5)

name="TEST_mix_pipi00"




G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="output/"


function sim_factory(modB)
    if modB < Bmin
        return SimulationParameters(name,
        lattice=geom.PyroPrimitive(1,1,1),
        A=[0 0 0 0],
        Jpm=Jpm,
        B=[1.,1.,0.]/√3,
        nsample=1000,
        kappa=2.0
        )

    else
        return SimulationParameters(name,
        lattice=geom.PyroPrimitive(2,1,1),
        A =  construct_landau_gauge(lat, [0 0 0 π; 0 0 0 0; 0 0 0 0]),
        Jpm=Jpm,
        B=[1.,1.,0.]/√3,
        nsample=1000,
        kappa=2.0
        )
    end
end



f = integrated_fieldsweep(data_dir, 
                          sim_factory=sim_factory,
                          magnetic_field_strengths=collect(range(0,0.5,50)),
                          ip=integration_settings["slow"],
                          Egrid=Egrid, g_tensor=G
                         )
println("Saving data to $(f)")





