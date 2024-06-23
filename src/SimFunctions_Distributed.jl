using Distributed

@everywhere begin
	include("./IOFunctions.jl")

end

@everywhere begin
	using SharedArrays
	using ProgressMeter

end



function calc_spinons_along_path(output_dir;
    sim::SimulationParameters,
    λ::Float64,
    path::BZPath)
    num_K = length(path.K)
	bands = SharedArray{Float64}(num_K, length(sim.lat.tetra_sites))
    
	
	status = @showprogress pmap(1:num_K) do J
        bands[J,:] = begin
            try
                spinon_dispersion(path.K[J], sim, λ)[1]'
            catch e
                if e isa DomainError
                    NaN*ones(Float64, 1,length(sim.lat.tetra_sites))
                else
                    throw(e)
                end
            end
        end
    end

    return save_spinons(output_dir;
		bands=sdata(bands),
        BZ_path=path,
        sim=sim)
end


"""
calc_spectral_weight_along_path(sim::SimulationParameters,
    ip::IntegrationParameters, Egrid::Vector{Float64}, path::BZPath)

    @param sim the simulation parameters
    @param ip the numerical constants, passed directly to spectral_weight
    @param Egrid the grid for energy values
    @param path the path through the Brillouin zone
    @param g_tensor the matrix such that, with respect to local axes, m = g S

    Saves the results to
"""
function calc_spectral_weight_along_path(
    output_dir::String;
    sim::SimulationParameters, λ::Float64,
    ip::IntegrationParameters, Egrid::Vector{Float64}, path::BZPath,
    g_tensor=nothing
    )
    
    num_K = length(path.K)
    
	Spm = SharedArray{ComplexF64}(num_K, length(Egrid))
	Spp = SharedArray{ComplexF64}(num_K, length(Egrid))
	Smagnetic = SharedArray{Float64}(num_K, length(Egrid))
	bounds = SharedArray{Float64}(num_K, 2)
    

    tmp_intensity = Sqω_set(Egrid)


	@showprogress @distributed for I=1:num_K
        k = path.K[I]*0.5
        q = SVector(k[1], k[2], k[3])

        try
            # no race condition here, yay
            spectral_weight!(tmp_intensity, q, sim, λ, ip, g_tensor )
            Spm[I, :] .= tmp_intensity.Sqω_pm
            Spp[I, :] .= tmp_intensity.Sqω_pp
            Smagnetic[I, :] .= tmp_intensity.Sqω_magnetic
            bounds[I,:] .= tmp_intensity.bounds
        catch e
            if e isa DomainError
                println("Q= $(q), negative dispersion: $(e)")
            else
                throw(e)
            end
        end
    end

    
    # save the data
    return save_SQW(output_dir, 
		Spm=sdata(Spm),
		Spp=sdata(Spp),
		Smagnetic=sdata(Smagnetic), 
		bounds=sdata(bounds),
        BZ_path=path,
        Egrid=collect(Egrid),
        sim=sim, 
        ip=ip
        )
end


"""
function calc_integrated_specweight(
    output_dir::String;
    sim::SimulationParameters, λ::Float64,
    ip::IntegrationParameters, Egrid::Vector{Float64},g_tensor=nothing)

Calculates the BZ integrated spectral weight for a particular sim, 
    saving results to output_dir
"""
function calc_integrated_specweight(
    output_dir::String;
    sim::SimulationParameters, λ::Float64,
    ip::IntegrationParameters, Egrid::Vector{Float64},g_tensor=nothing)


#=
    Sω_pm = zeros(ComplexF64,size(Egrid))
    Sω_pp = zeros(ComplexF64,size(Egrid))
    Sω_magnetic = zeros(Float64,size(Egrid))
=#



	res = @showprogress @distributed ((x,y)->x.+y) for _ in 1:ip.n_K_samples
            q = (1 .- 2 .*(@SVector rand(3)))*4π/8
            p = (1 .- 2 .*(@SVector rand(3)))*4π/8

            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, λ, g_tensor)

            (
			broadened_peaks( S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE ),
            broadened_peaks( S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE ),
            broadened_peaks( S_magnetic_rs::Matrix{Float64}, E_rs, Egrid, ip.broadening_dE )
			)
    end



    


    (Spm_res, Spp_res, Smagnetic_res) = res

    # add the threads' results together and save the data
    return save_SQW(output_dir, 
                Spm=Spm_res,
                Spp=Spp_res,
                Smagnetic=Smagnetic_res,
                bounds=nothing, BZ_path=nothing, Egrid=collect(Egrid),
                sim=sim, ip=ip, prefix="integrated")
end




