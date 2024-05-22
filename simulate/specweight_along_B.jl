include("../src/SimFunctions.jl")

magnetic_fields = map( b-> b*[1,1,1]/√3, range(0,0.3, 50) )

simlist = map(
    b->SimulationParameters("piflux-stretch",
    A=[ 0 0 π π ; 0 0 0 0; 0 0 π π; 0 0 0 0 ],
    Jpm=1/3,
    B=b,
    nsample=1000,
    kappa=2.0
    ),
    magnetic_fields);


G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="output_b_vary/"

for sim in simlist
    f = calc_integrated_specweight(data_dir, 
	    sim=sim,
        ip=integration_settings["fast"],
	    Egrid=Egrid, g_tensor=G
       )
    println("Saving data to $(f)")
end
