using Distributed

include("driver-headless.jl")

magnetic_fields = [
    [0.,0.,0.], [1,1,1]*0.06/√3, [0,1,1]*0.06/√2 
    ]

simlist_X = map(
    b->SimulationParameters("TEST_0flux",
    A=zeros(Float64, 4,4),
    Jpm=-0.046,
    B=b,
    nsample=1000,
    kappa=2.0
    ),
    magnetic_fields);

simlist = [SimulationParameters(s, 0.) for s in simlist_X]

for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output/",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["very_fast"]
        )
end


    
