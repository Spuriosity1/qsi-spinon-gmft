include("driver-headless.jl")

magnetic_fields = [
    [0.,0.,0.], [1,1,1]*0.1/√3, [0,1,1]*0.1/√2 
    ]

simlist = map(
    b->SimulationParameters("piflux-2",
    A=[ 0 0 π π ; 0 0 0 0; 0 0 π π; 0 0 0 0 ],
    Jpm=1/3,
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


    
