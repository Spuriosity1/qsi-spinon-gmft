using Distributed

include("driver.jl")


Jpm = -0.05

Bmin =sqrt(-9*Jpm/5)
println("minimum B = ",Bmin)

magnetic_fields = [
bb*[0,1,1]/√3 for bb in [Bmin, Bmin*1.5, Bmin*2]
    ]

simlist_X = map(
    b->SimulationParameters("piflux+1e-4",
    A =  [ 0 0 0 π ; 0 0 0 0; 0 0 0 π; 0 0 0 0 ],
    Jpm=Jpm,
    B=b,
    nsample=1000,
    kappa=2.0
    ),
    magnetic_fields);

simlist = [ SimulationParameters(s, 1e-4) for s in simlist_X ]

for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output/",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["very_slow"]
        )
end


    
