include("../src/SimFunctions.jl")

Egrid=collect(range(0,3.2,150))


Jpm = -0.05

Bmin =sqrt(-9*Jpm/5)

name="TEST_mix_pipi00"

sim_0flux = SimulationParameters(name,
    A =  [ 0 0 0 0 ; 0 0 0 0; 0 0 0 0; 0 0 0 0 ],
    Jpm=Jpm,
    B=[1.,1.,0.]/√3,
    nsample=1000,
    kappa=2.0
    )



sim_ππ00 = SimulationParameters(name,
    A =  [ 0 0 0 π ; 0 0 0 0; 0 0 0 π; 0 0 0 0 ],
    Jpm=Jpm,
    B=[1.,1.,0.]/√3,
    nsample=1000,
    kappa=2.0
    )


G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

data_dir="output/"


function sim_factory(modB)
    if modB < Bmin
        sim = sim_0flux
    else
        sim = sim_ππ00
    end
    return SimulationParameters(sim.name,
                                        A=sim.A, Jpm=sim.Jpm,
                                        B=modB*[1.,1.,0.]/√3,
                                        nsample = 10000,
                                        kappa=2.0)
end

f = integrated_fieldsweep(data_dir, 
                          sim_factory=sim_factory,
                          magnetic_field_strengths=collect(range(0,0.5,50)),
                          ip=integration_settings["very_fast"],
                          Egrid=Egrid, g_tensor=G
                         )
println("Saving data to $(f)")





