include("../src/SimFunctions.jl")


sim = SimulationParameters("piflux-along",
    A=[ 0 0 π π ; 0 0 0 0; 0 0 π π; 0 0 0 0 ],
    Jpm=1/3,
    B=[1.,1.,1.]/√3,
    nsample=1000,
    kappa=2.0
    )


G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="output/"

f = integrated_fieldsweep(data_dir, 
                          sim=sim,
                          magnetic_field_strengths=collect(range(0,0.1,50)),
                          ip=integration_settings["fast"],
                          Egrid=Egrid, g_tensor=G
                         )
println("Saving data to $(f)")





