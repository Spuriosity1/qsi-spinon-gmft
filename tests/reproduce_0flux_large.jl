include("../simulate/driver-headless.jl")
include("../src/PlotFunctions.jl")


sim = SimulationParameters("0flux",
    lattice=geom.PyroPrimitive(2,2,1),
    A=[0 0 0 0; 0 0 0 0; 0 0 0 0; 0 0 0 0],
    Jpm=-0.046,
    B=[0.,0.,0.],
    n_samples=10000
    )


df = run_sim(
        data_dir="output",
        figure_dir="",
        sim=sim, 
        integral_params=integration_settings["fast"],
        calc_integrated=false
        )


print(df)

plot_spinons(load(df["spinon"]))
savefig("output/$(sim.name)-spinon.png" )
plot_spectral_weight(load(df["specweight"]))
savefig("output/$(sim.name)-specweight.png" )
