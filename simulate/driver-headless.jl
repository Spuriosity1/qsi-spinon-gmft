include("../src/SimFunctions.jl")

function run_sim(;data_dir, figure_dir,
	sim::SimulationParameters,
	integral_params::IntegrationParameters,
	calc_specweight=true,
	calc_integrated=true,
	k_density_spinon_dispersion=60,
	k_density_specweight=20
	)


	G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

	path_spinons = generate_path(geom.high_symmetry_points, 
	    split("\\Gamma X W K \\Gamma L U W"), 
		points_per_unit=k_density_spinon_dispersion, K_units=2π/8)
	
	path = generate_path(geom.high_symmetry_points, 
	    split("\\Gamma X W K \\Gamma L U W"),
		points_per_unit=k_density_specweight, K_units=4π/8)
	
	println("Computing spinon dispersions...")
	# compute spinons
	d = calc_spinons_along_path(data_dir, sim=sim, path=path_spinons)
		
	println("Calculating spectral weight data...")
	datafiles = []

	# autorange this based on the spinon dispersion
	max_E = 2.2* maximum( load(d)["spinon_dispersion"]["bands"] )
	Egrid = collect(range(0,max_E,150)) # TODO consider updating this based on broadening_dE
	
	
	# k-resolved simulation
	if calc_specweight
	    @printf("Running specweight simulation...\n")
	    f = calc_spectral_weight_along_path(data_dir, 
	    sim=sim,
	    ip=integral_params, 
	    Egrid=Egrid, path=path, g_tensor=G)
	    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
	    push!(datafiles, f)
	    println("Saving data to ",f)
	end
	
		
	# BZ integrated
	if calc_integrated
	    @printf("Running BZ integrated specweight simulation...\n")
	    f = calc_integrated_specweight(data_dir, 
	    sim=sim,
	    ip=integral_params, 
	    Egrid=Egrid, g_tensor=G)
	    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
	    push!(datafiles, f)
	    println("Saving data to ",f)
	end
	
end		
