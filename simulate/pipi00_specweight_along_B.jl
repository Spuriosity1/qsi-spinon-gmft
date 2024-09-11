include("../src/SimFunctions.jl")

Egrid=collect(range(0,3.2,150))


Jpm = -0.04

Bmin =sqrt(-9*Jpm/5)

name="mix_pipi00"


G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="../output/"

field_direction = (@SVector [1.,1.,0.])/√3


function sim_factory(modB)
    if modB < Bmin
        return SimulationParameters(name,
        lattice=geom.PyroPrimitive(1,1,1),
        A=[0 0 0 0],
        Jpm=Jpm,
        B=modB*[1.,1.,0.]/√3,
        nsample=1000,
        kappa=2.0
        )

    else
        return SimulationParameters(name,
        lattice=geom.PyroPrimitive(2,1,1),
        A =  construct_landau_gauge(lat, [0 0 0 π; 0 0 0 0; 0 0 0 0]),
        Jpm=Jpm,
        B=modB*[1.,1.,0.]/√3,
        nsample=1000,
        kappa=2.0
        )
    end
end



ip = IntegrationParameters(n_K_samples=100000,broadening_dE=0.02)
# integration_settings["very_slow"]

f = integrated_fieldsweep(data_dir, 
                          sim_factory=sim_factory,
                          magnetic_field_strengths=collect(range(0,0.5,50)),
                          ip=ip,
                          Egrid=Egrid, g_tensor=G
                         )
println("Saving data to $(f)")





