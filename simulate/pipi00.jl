include("../src/PlotFunctions.jl")

A_ππ00 = [ 0 0 0 π ; 0 0 0 0; 0 0 0 π; 0 0 0 0 ]


Jpm = -0.05

# phase transition 0-> ππ00 is at B = √(9/5 * -Jpm) 
Bmin =sqrt(-9*Jpm/5)
println("minimum B = ",Bmin)

# a cheap hack to see different directions
magnetic_fields = [B0*[1,1,0]/sqrt(2)  for B0 in (Bmin, Bmin+0.1, Bmin+0.2,Bmin+0.3)]
    

println("Calculating chemical potentials...")

simlist_x = []
@showprogress for b in magnetic_fields
    push!(simlist_x,
        SimulationParameters("pipi00-kludge", A=A_ππ00, Jpm=Jpm, B=b, nsample=10000, kappa=2.0)
        )
end

simlist = [SimulationParameters(s, 1e-3) for s in simlist_x]


ip = integration_settings["very_slow"]

Egrid = collect(range(0,1.4,150))

G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

figure_dir = "figures/"
data_dir = "output/"



path_spinons = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=60, K_units=2π/8)

println("Plotting spinon dispersions...")
# plot the spinons
@showprogress for sim in simlist
	d = calc_spinons_along_path(data_dir, sim=sim, path=path_spinons)
	p = plot_spinons(load(d))
    savefig(p, figure_dir*"spinon_dispersion"*sim_identifier(sim)*".pdf")
end

path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=30, K_units=4π/8)

println("Calculating spectral weight data...")
datafiles = []
# run the simulation

for (j,sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", j, length(simlist))
    f = calc_spectral_weight_along_path(data_dir, 
    sim=sim,
    ip=ip, 
    Egrid=Egrid, path=path, gtensor=G)
    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
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

