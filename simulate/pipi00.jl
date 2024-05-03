include("../src/PlotFunctions.jl")

A_ππ00 = [ 0 0 0 π ; 0 0 0 0; 0 0 0 π; 0 0 0 0 ]


Jpm = -0.05

# phase transition 0-> ππ00 is at B = √(9/5 * -Jpm) 
Bmin =sqrt(-9*Jpm/5)
println("minimum B = ",Bmin)

simlist = [SimulationParameters("pipi00", A=A_ππ00, Jpm=Jpm, B=B0*[1,1,0]/sqrt(2), nsample=1000, kappa=2.0) for B0 in (Bmin, Bmin+0.1, Bmin+0.2,Bmin+0.3)] 

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

println("Calculating the averaged spectral weight")


