
using Pkg 
Pkg.activate(joinpath(@__DIR__,"../"))
Pkg.instantiate()
Pkg.precompile()


using LinearAlgebra
using StaticArrays
using ProgressMeter
using JLD
using HDF5
using Printf


include("BZMath.jl")
include("SpinonStructure.jl")

using .BZmath
using .SpinonStructure


    
function save_SQW(output_dir::String;
        Spm::Union{Matrix{ComplexF64},Vector{ComplexF64}},
        Spp::Union{Matrix{ComplexF64},Vector{ComplexF64}},
        Smagnetic::Union{Matrix{Float64},Vector{Float64}},
        bounds::Union{Nothing,Matrix{Float64}},
        Egrid::Vector{Float64},
        BZ_path::Union{Nothing,BZPath},
        sim::SimulationParameters, 
        ip::IntegrationParameters,
        prefix="SQW")

    name = output_dir*"/"*prefix*sim_identifier(sim)*".jld"
    rm(name, force=true)
    jldopen(name, "w") do file
        g1 = create_group(file, "integration_parameters")
        g1["n_K_samples"] = ip.n_K_samples
        g1["broadening_dE"] = ip.broadening_dE
        g1["version"] = 1.0
#        g1["BZ_grid_density"] = ip.BZ_grid_density

        g = create_group(file, "physical_parameters")
        g["name"] = sim.name
        g["emergent_fluxes"] = calc_fluxes(sim.A)
        g["gauge"] = sim.A
        g["Jpm"] = sim.Jpm
        g["B"] = sim.B
        g["lambda"] = sim.λ
        g["L"]=sim.lat.L

        d = create_group(file, "intensity") 
        d["Spm"] = Spm
        d["Spp"] = Spp
        d["Smagnetic"] = Smagnetic
        if bounds != nothing
            d["bounds"] = bounds
        end
        # a list of K points, such that the I'th S slics corresponds to the I'th K point
        if BZ_path != nothing
            d["Q_list"] = BZ_path.K 
            d["tau"] = BZ_path.t
            d["ticks_tau"] = BZ_path.ticks_t
            d["ticks_label"] = BZ_path.ticks_label
        end
        d["W"] = Egrid
    end
    return name
end


function save_field_sweep(output_dir::String;
        Spm::Union{Matrix{ComplexF64},Vector{ComplexF64}},
        Spp::Union{Matrix{ComplexF64},Vector{ComplexF64}},
        Smagnetic::Union{Matrix{Float64},Vector{Float64}},
        Egrid::Vector{Float64},
        magnetic_field_strengths::Vector{Float64},
        sim::SimulationParameters, 
        ip::IntegrationParameters,
        prefix="int_fieldsweep")

    name = output_dir*"/"*prefix*sim_identifier(sim)*".jld"
    rm(name, force=true)
    jldopen(name, "w") do file
        g1 = create_group(file, "integration_parameters")
        g1["n_K_samples"] = ip.n_K_samples
        g1["broadening_dE"] = ip.broadening_dE
        g1["version"] = 1.0
#        g1["BZ_grid_density"] = ip.BZ_grid_density

        g = create_group(file, "physical_parameters")
        g["name"] = sim.name
        g["emergent_fluxes"] = calc_fluxes(sim.A)
        g["gauge"] = sim.A
        g["Jpm"] = sim.Jpm
        g["B"] = sim.B
        g["lambda"] = sim.λ
        g["L"]=sim.lat.L

        d = create_group(file, "integrated_intensity") 
        d["Spm"] = Spm
        d["Spp"] = Spp
        d["Smagnetic"] = Smagnetic
        d["W"] = Egrid
        d["magnetic_field_strengths"] = magnetic_field_strengths
    end
    return name
end



function save_spinons(output_dir::String;
        bands::Matrix{Float64},
        BZ_path::BZPath,
        sim::SimulationParameters, 
        prefix="spinons")

    name = output_dir*"/"*prefix*sim_identifier(sim)*".jld"
    jldopen(name, "w") do file

        g = create_group(file, "physical_parameters")
        g["name"] = sim.name
        g["fluxes"] = sim.A
        g["Jpm"] = sim.Jpm
        g["B"] = sim.B
        g["lambda"] = sim.λ
        g["L"]=sim.lat.L

        d = create_group(file, "spinon_dispersion") 
        d["bands"] = bands
        # a list of K points, such that the I'th S slics corresponds to the I'th K point
        d["Q_list"] = BZ_path.K 
        d["tau"] = BZ_path.t
        d["ticks_tau"] = BZ_path.ticks_t
        d["ticks_label"] = BZ_path.ticks_label
    end
    return name
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
    sim::SimulationParameters,
    ip::IntegrationParameters, Egrid::Vector{Float64}, path::BZPath,
    g_tensor=nothing
    )
    
    num_K = length(path.K)
    
    Spm = zeros(ComplexF64, num_K, length(Egrid))
    Spp = zeros(ComplexF64, num_K, length(Egrid))
    Smagnetic = zeros(Float64, num_K, length(Egrid))
    bounds = zeros(Float64, num_K, 2)
    
    p = Progress(num_K,desc="Spectral weight: ")

    Threads.@threads for I = 1:num_K
        k = path.K[I]*0.5
        q = SVector(k[1], k[2], k[3])
        # hard-coded DO g-tensor
        try
            # no race condition here, yay
            Spm[I, :], Spp[I, :], Smagnetic[I, :], bounds[I,:] = spectral_weight(q, Egrid, sim, ip, g_tensor )
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



function calc_integrated_specweight(
    output_dir::String;
    sim::SimulationParameters,
    ip::IntegrationParameters, Egrid::Vector{Float64},g_tensor=nothing)


#=
    Sω_pm = zeros(ComplexF64,size(Egrid))
    Sω_pp = zeros(ComplexF64,size(Egrid))
    Sω_magnetic = zeros(Float64,size(Egrid))
=#

    prog = Progress(ip.n_K_samples, desc="BZ Average")
 
    Spm_res = zeros(ComplexF64, Threads.nthreads(), length(Egrid) ) 
    Spp_res = zeros(ComplexF64, Threads.nthreads(), length(Egrid) ) 
    Smagnetic_res = zeros(Float64, Threads.nthreads(), length(Egrid) ) 

    Threads.@threads for _ = 1:ip.n_K_samples  
        q = (1 .- 2 .*(@SVector rand(3)))*4π/8
        p = (1 .- 2 .*(@SVector rand(3)))*4π/8

        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, g_tensor)
            
            # race condition time
            id = Threads.threadid()
            Spm_res[id, :] += broadened_peaks(S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            Spp_res[id, :] += broadened_peaks(S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            Smagnetic_res[id, :] += broadened_peaks(S_magnetic_rs::Matrix{Float64}, E_rs, Egrid,
                                         ip.broadening_dE )
                    
        catch e
            if e isa DomainError
                println("p=$(p) or p-q=$(p-q): Negative dispersion")
            else
                throw(e)
            end
        end
        next!(prog)
    end
    

    finish!(prog)


    # add the threads' results together and save the data
    return save_SQW(output_dir, 
                    Spm=sum(eachrow(Spm_res)),
                    Spp=sum(eachrow(Spp_res)),
                    Smagnetic=sum(eachrow(Smagnetic_res)),
                    bounds=nothing,
                    BZ_path=nothing,
                    Egrid=collect(Egrid), sim=sim, ip=ip,
                    prefix="integrated"
                   )
end





# single-threaded version
function calc_integrated_S(
    output_dir::String;
    sim::SimulationParameters,
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
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, g_tensor) 
            Spm_res += broadened_peaks(S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            Spp_res += broadened_peaks(S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, ip.broadening_dE )
            Smagnetic_res += broadened_peaks(S_magnetic_rs::Matrix{Float64}, E_rs, Egrid,
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
    sim_factory, # A map taking in a bield strength and returning SimulationParameters
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

        Spm_res, Spp_res, Smagnetic_res = calc_integrated_S(output_dir,
                                                            sim=this_sim,
                                                            ip=ip,
                                                            Egrid=Egrid,
                                                            g_tensor=g_tensor
                                                           )
        Spm[J,:] = Spm_res
        Spp[J,:] = Spp_res
        Smagnetic[J,:] = Smagnetic_res

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




function calc_spinons_along_path(output_dir;
    sim::SimulationParameters,
    path::BZPath)
    num_K = length(path.K)
    bands = zeros(Float64, num_K, length(sim.lat.tetra_sites))
    
    prog=Progress(num_K, desc="Spinon Dispersion: ") 
    @Threads.threads for J=1:num_K
        bands[J,:] = begin
            try
                spinon_dispersion(path.K[J], sim )[1]'
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


# some useful defaults

const integration_settings = Dict(
    "very_fast" => IntegrationParameters(n_K_samples=10,   broadening_dE=0.1),
    "fast" =>   IntegrationParameters(n_K_samples=100,  broadening_dE=0.05),
    "slow" =>      IntegrationParameters(n_K_samples=1000, broadening_dE=0.02),
    "very_slow" => IntegrationParameters(n_K_samples=10000,broadening_dE=0.02)
)

