include("../src/PlotFunctions.jl")

A = zeros(Float64, 4,4)

magnetic_fields = vcat(
    map(b->(b/√3)*[1.,1.,1.], [0,0.05]),
)


simlist = map(
    b->SimulationParameters("TEST_0flux",
        A=A,
        Jpm=-0.046,
        B=b,
        nsample=1000,
        kappa=2.0
        ),
    magnetic_fields);

ip = integration_settings["fast"]

Egrid = collect(range(0,1.4,150))

G = [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]


figure_dir = "figures/"
data_dir = "output/"

path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=30, K_units=4π/8)

println("Plotting spinon dispersions...")
# plot the spinons
@showprogress for sim in simlist
    p = plot_spinons(sim,path)
    savefig(p, figure_dir*"spinon_dispersion"*sim_identifier(sim)*".pdf")
end


println("Calculating spectral weight data...")
datafiles = []
# run the simulation
for (j,sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", j, length(simlist))
    f = calc_spectral_weight_along_path(sim, ip, Egrid, path, G, data_dir)
    push!(datafiles, f)
    println("Saving data to ",f)
end

println("Plotting the spectral weights")
for specweight_data in datafiles
    data = load(specweight_data)
    sim = SimulationParameters(data["physical_parameters"])

    p = plot_spectral_weight(data,"Spm")
    savefig(p, figure_dir*"corr_S+-"*sim_identifier(sim)*".pdf")


    p = plot_spectral_weight(data,"Spp")
    savefig(p, figure_dir*"corr_S++"*sim_identifier(sim)*".pdf")


    p = plot_spectral_weight(data,"Smagnetic")
    savefig(p, figure_dir*"corr_Smagnetic"*sim_identifier(sim)*".pdf")
end

