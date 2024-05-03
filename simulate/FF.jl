include("../src/PlotFunctions.jl")

A_FF =  load_A("gaugefiles/FF_even_pi18.gauge")

"""
Returns the ratio of g0 to g1 in the special case B || [111], given a particular value of Φ0
"""
g0_g123_ratio_from_flux(desired_Φ0) = (1-4*cos(desired_Φ0/3)^2)^-1

"""
FF_111_B(desired_Φ, Jpm)
-> B

x = g0/g1 <0

Returns the (111) field needed to realise the prescribed field
"""
function FF_111_B(desired_Φ0, Jpm)
    # x = Φ0 -> ( 1 - 2*cos((π/2 - Φ0)*2/3) )^-1;
    xp = g0_g123_ratio_from_flux(desired_Φ0)
    
    B = sqrt(Jpm* 6*(1-xp)/ (xp*5/9 - 5) )
    return B
end


Φ0 = 2π/3

simlist = map(
    Jpm->SimulationParameters("FF",
        A=A_FF, 
        Jpm=Jpm,
        B=FF_111_B(Φ0, Jpm)* [1.,1.,1.]/sqrt(3),
        nsample=1000,
        kappa=2.0
    ),
    [-0.01, -0.02, -0.05, -0.1]
)


ip = integration_settings["very_fast"]

Egrid = collect(range(0,3,150))

figure_dir = "figures/"

path = generate_path(geom.high_symmetry_points, 
    split("\\Gamma X W K \\Gamma L U W"), points_per_unit=5, K_units=4π/8)

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

