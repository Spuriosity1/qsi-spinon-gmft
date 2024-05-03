include("../src/PlotFunctions.jl")

A =  [ 0 0 π π ; 0 0 0 0; 0 0 π π; 0 0 0 0 ]


magnetic_fields = vcat(
    map(b->(b/√3)*[1.,1.,1.], [0,0.01,0.02,0.03,0.04,0.05]),
    map(b->(b/√2)*[1.,1.,0.], [0,0.01,0.02,0.03,0.04,0.05])
)

simlist = map(
    b->SimulationParameters("piflux",
        A=A,
        Jpm=1./3,
        B=b,
        nsample=1000,
        kappa=2.0
        ),
    magnetic_fields);

ip = integration_settings["very_slow"]

Egrid = collect(range(0,3,150))

figure_dir = "figures/"

path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=10, K_units=4π/8)

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
    f = calc_spectral_weight_along_path(sim, ip, Egrid, path, figure_dir)
    push!(datafiles, f)
end

println("Plotting the spectral weights")
for specweight_data in datafiles
    data = load(specweight_data)
    p = plot_spectral_weight(data)
    sim = SimulationParameters(data["physical_parameters"])
    savefig(p, figure_dir*"spectral_weight"*sim_identifier(sim)*".pdf")
end

