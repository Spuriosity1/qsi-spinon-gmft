
include("driver-headless.jl")

#magnetic_fields = [
#    [0.,0.,0.], [1,1,1]*0.06/√3, [0,1,1]*0.06/√2 
#    ]


magnetic_fields = [
    [1,1,1]*0.1/√3, [0,1,1]*0.1/√2,
    [1,1,1]*0.2/√3, [0,1,1]*0.2/√2,
    [1,1,1]*0.3/√3, [0,1,1]*0.3/√2 
    ]

simlist = map(
    b->SimulationParameters("0flux",
    A=zeros(Float64, 4,4),
    Jpm=-0.046,
    B=b,
    nsample=1000,
    kappa=2.0
    ),
    magnetic_fields);


for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output/",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["very_slow"]
        )
end


    