include("../src/PlotFunctions.jl")

A_ππ00 = [ 0 0 0 π ; 0 0 0 0; 0 0 0 π; 0 0 0 0 ]


Jpm = -0.05

# phase transition 0-> ππ00 is at B = √(9/5 * -Jpm) 
Bmin =sqrt(-9*Jpm/5)
println("minimum B = ",Bmin)

# a cheap hack to see different directions
magnetic_fields = [B0*[1,1,0]/sqrt(2)  for B0 in (Bmin, Bmin+0.1, Bmin+0.2,Bmin+0.3)]
    

println("Calculating chemical potentials...")

simlist = []
@showprogress for b in magnetic_fields
    push!(simlist,
        SimulationParameters("pipi00", A=A_ππ00, Jpm=Jpm, B=b, nsample=1000000, kappa=2.0)
        )
end

ip = integration_settings["very_slow"]

Egrid = collect(range(0,3,150))

G = [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

figure_dir = "figures/"
data_dir = "output/"



path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=60, K_units=4π/8)

println("Plotting spinon dispersions...")
# plot the spinons
@showprogress for sim in simlist
    p = plot_spinons(sim,path)
    savefig(p, figure_dir*"spinon_dispersion"*sim_identifier(sim)*".pdf")
end


path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=30, K_units=4π/8)

println("Calculating spectral weight data...")
datafiles = []
# run the simulation
for (j,sim) in enumerate(simlist)
    #f = calc_spectral_weight_along_path(sim, ip, Egrid, path, G, data_dir)

    f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
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

