include("../src/PlotFunctions.jl")

A = zeros(Float64, 4,4)

magnetic_fields = [
    [0.,0.,0.], [1,1,1]*0.06/√3, [0,1,1]*0.06/√2 
    ]

simlist = map(
    b->SimulationParameters("0flux",
        A=A,
        Jpm=-0.046,
        B=b,
        nsample=1000,
        kappa=2.0
        ),
    magnetic_fields);

ip = integration_settings["very_slow"]

Egrid = collect(range(0,1.4,150))

G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

figure_dir = "figures/"
data_dir = "output/"


path_spinons = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=60, K_units=2π/8)

path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=30, K_units=4π/8)

println("Plotting spinon dispersions...")
# plot the spinons
@showprogress for sim in simlist
	d = calc_spinons_along_path(data_dir, sim=sim, path=path_spinons)
	p = plot_spinons(load(d))
    savefig(p, figure_dir*"spinon_dispersion"*sim_identifier(sim)*".pdf")
end


println("Calculating spectral weight data...")
datafiles = []


# k-resolved simulation
#=
for (j,sim) in enumerate(simlist)
    @printf("Running specweight simulation %d of %d\n", j, length(simlist))
    f = calc_spectral_weight_along_path(data_dir, 
    sim=sim,
    ip=ip, 
    Egrid=Egrid, path=path, g_tensor=G)
    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
    push!(datafiles, f)
    println("Saving data to ",f)
end
=#

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

# BZ integrated
for (j,sim) in enumerate(simlist)
    @printf("Running powder simulation %d of %d\n", j, length(simlist))
    f = calc_integrated_specweight(data_dir, 
    sim=sim,
    ip=ip, 
    Egrid=Egrid, g_tensor=G)
    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
    push!(datafiles, f)
    println("Saving data to ",f)
end
