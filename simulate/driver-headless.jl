using Pkg; Pkg.activate(".")
include("../src/SimFunctions.jl")

function run_sim(;data_dir, figure_dir,
	sim::SimulationParameters,
	integral_params::IntegrationParameters,
	calc_specweight=true,
	calc_integrated=true,
	k_density_spinon_dispersion=60,
	k_density_specweight=20,
	fudge_lambda=0.
	)


	G = @SMatrix [0. 0. 0.;
	 0. 0. 0.;
	 1. 0. 0.]

	 dfiles = Dict()

	path = generate_path(geom.high_symmetry_points, 
	    split("\\Gamma X W K \\Gamma L U W"), 
		points_per_unit=k_density_spinon_dispersion, K_units=2π/8)


    println("Calculating large-N spinon mass")
	csim = CompiledModel(sim)
    
	println("Computing spinon dispersions...")
	# compute spinons
	d = calc_spinons_along_path(data_dir, csim=csim, path=path)
		
	datafiles = []

	# autorange this based on the spinon dispersion
	band_data = load(d)["spinon_dispersion"]["bands"]
	band_data[ isnan.(band_data) ] .= - Inf # an ugly hack
	max_E = 2.2* maximum( band_data )
	Egrid = collect(range(0,max_E,150)) # TODO consider updating this based on broadening_dE
    println("Max energy for specweight: $(max_E)")
	dfiles["spinon"] = d
	
	
	# k-resolved simulation
	if calc_specweight
	    @printf("Running specweight simulation...\n")
	    f = calc_spectral_weight_along_path(data_dir, 
	    sim=csim,
	    ip=integral_params, 
	    Egrid=Egrid, path=path, g_tensor=G)
	    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
	    push!(datafiles, f)
	    println("Saving data to ",f)

		dfiles["specweight"] = f
	end
	
		
	# BZ integrated
	if calc_integrated
	    @printf("Running BZ integrated specweight simulation...\n")
	    f = calc_integrated_specweight(data_dir, 
	    csim=csim,
	    ip=integral_params, 
	    Egrid=Egrid, g_tensor=G)
	    # f = data_dir*"/SQW"*sim_identifier(sim)*".jld"
	    push!(datafiles, f)
	    println("Saving data to ",f)


		dfiles["integrated"] = f
	end

	return dfiles
	
end		
