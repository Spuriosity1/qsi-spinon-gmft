include("./IOFunctions.jl")


using ProgressMeter

function calc_spinons_along_path(output_dir;
    sim::SimulationParameters,
    λ::Float64,
    path::BZPath)
    num_K = length(path.K)
    bands = zeros(Float64, num_K, length(sim.lat.tetra_sites))
    
    prog=Progress(num_K, desc="Spinon Dispersion: ") 
    @Threads.threads for J=1:num_K
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
        next!(prog)
    end
    finish!(prog)

    return save_spinons(output_dir;
        bands=bands,
        BZ_path=path,
        sim=sim)
end

###########
# expensive calculations
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
    
    Spm = zeros(ComplexF64, num_K, length(Egrid))
    Spp = zeros(ComplexF64, num_K, length(Egrid))
    Smagnetic = zeros(Float64, num_K, length(Egrid))
    bounds = zeros(Float64, num_K, 2)
    
    p = Progress(num_K,desc="Spectral weight: ")

    tmp_intensity = Sqω_set(Egrid)


    Threads.@threads for I = 1:num_K
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
        next!(p)
    end
    finish!(p)

    
    # save the data
    return save_SQW(output_dir, 
        Spm=Spm,
        Spp=Spp,
        Smagnetic=Smagnetic, 
        bounds=bounds,
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

    prog = Progress(ip.n_K_samples, desc="BZ Average")

    samples_per_trial = ip.n_K_samples ÷ Threads.nthreads()

    
    function single_trial()        
        Spm_res = zeros(ComplexF64, length(Egrid) ) 
        Spp_res = zeros(ComplexF64, length(Egrid) ) 
        Smagnetic_res = zeros(Float64, length(Egrid) ) 


        for _ in 1:samples_per_trial
            q = (1 .- 2 .*(@SVector rand(3)))*4π/8
            p = (1 .- 2 .*(@SVector rand(3)))*4π/8

            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, λ, g_tensor)

            broadened_peaks!(Spm_res, S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            broadened_peaks!(Spp_res, S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            broadened_peaks!(Smagnetic_res, S_magnetic_rs::Matrix{Float64}, E_rs, Egrid, ip.broadening_dE )

            next!(prog)
        end
        return (Spm_res, Spp_res, Smagnetic_res)
    end

    tasks = map(1:Threads.nthreads()) do chunk
        Threads.@spawn single_trial()
    end

    chunk_sums = fetch.(tasks)
    
    finish!(prog)


    (Spm_res, Spp_res, Smagnetic_res) = reduce((x,y)->x.+y, chunk_sums)

    # add the threads' results together and save the data
    return save_SQW(output_dir, 
                Spm=Spm_res,
                Spp=Spp_res,
                Smagnetic=Smagnetic_res,
                bounds=nothing, BZ_path=nothing, Egrid=collect(Egrid),
                sim=sim, ip=ip, prefix="integrated")
end





# single-threaded version
function calc_integrated_S(
    output_dir::String;
    sim::SimulationParameters,
    λ::Float64,
    ip::IntegrationParameters,
    Egrid::Vector{Float64},
    g_tensor=nothing
    )
 
    Spm_res = zeros(ComplexF64, length(Egrid) ) 
    Spp_res = zeros(ComplexF64, length(Egrid) ) 
    Smagnetic_res = zeros(Float64, length(Egrid) ) 

    for _ = 1:ip.n_K_samples  
        q = (1 .- 2 .*(@SVector rand(3)))*4π/8
        p = (1 .- 2 .*(@SVector rand(3)))*4π/8

        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, λ, g_tensor) 
            Spm_res .+= broadened_peaks(S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            Spp_res .+= broadened_peaks(S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            Smagnetic_res .+= broadened_peaks(S_magnetic_rs::Matrix{Float64}, E_rs, Egrid,
                                         ip.broadening_dE )
                    
        catch e
            if e isa DomainError
                println("p=$(p) or p-q=$(p-q): Negative dispersion")
            else
                throw(e)
            end
        end
    end
    

    return Spm_res, Spp_res, Smagnetic_res
end

function integrated_fieldsweep(output_dir::String;
    sim_factory, # A map taking in a field strength and returning SimulationParameters
    magnetic_field_strengths::Vector{Float64},
    ip::IntegrationParameters,
    Egrid::Vector{Float64},
    g_tensor)

    num_B = length(magnetic_field_strengths)
    prog=Progress(num_B, desc="Field Sweep: ") 

    Spm = zeros(ComplexF64, num_B, length(Egrid))
    Spp = zeros(ComplexF64, num_B, length(Egrid))
    Smagnetic = zeros(Float64, num_B, length(Egrid))

    
    @Threads.threads for J=1:num_B
        this_sim = sim_factory(magnetic_field_strengths[J])::SimulationParameters
        this_λ = calc_lambda(this_sim)
        Spm_res, Spp_res, Smagnetic_res = calc_integrated_S(output_dir,
                                                            sim=this_sim,
                                                            λ=this_λ,
                                                            ip=ip,
                                                            Egrid=Egrid,
                                                            g_tensor=g_tensor
                                                           )

        Spm[J,:] .= Spm_res::Vector{ComplexF64}
        Spp[J,:] .= Spp_res::Vector{ComplexF64}
        Smagnetic[J,:] .= Smagnetic_res::Vector{Float64}

        next!(prog)
    end
    finish!(prog)


    save_field_sweep(output_dir::String;
        Spm=Spm,
        Spp=Spp,
        Smagnetic=Smagnetic,
        Egrid=Egrid,
        magnetic_field_strengths=magnetic_field_strengths,
        sim=sim_factory(magnetic_field_strengths[1]), 
        ip=ip
       )
end




