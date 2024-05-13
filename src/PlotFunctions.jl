using Distributed

@everywhere begin
  # instantiate environment
  using Pkg 
  Pkg.activate(joinpath(@__DIR__,"../"))
  Pkg.instantiate()
  Pkg.precompile()
end

@everywhere begin
using LinearAlgebra
using StaticArrays
using ProgressMeter
using Plots
using JLD
using HDF5
using Printf


include("BZMath.jl")
include("SpinonStructure.jl")
end

@everywhere begin
    using .BZmath
    using .SpinonStructure
end
#=
const kspace_points = Dict(
    "\\Gamma"=> [0.,0.,0.],
    "X"=> [1.,0.,0.],
    "Y"=> [0.,1.,0.],
    "Z"=> [0.,0.,1.],
    "W_X"=> [1.,0.5,0.],
    "W_Y"=> [0.,1.,0.5],
    "W_Z"=> [0.5,0.,1.],
    "K_{XY}"=> [0.75,0.75,0.],
    "K_{YZ}"=> [0.,0.75,0.75],
    "K_{ZX}"=> [0.75,0.,0.75],
    "L"=> [0.5,0.5,0.5],
    "U"=> [1.0, 0.25,0.25]
)
=#

# A 'hash function' providing an informative filenmae for simulation data
function sim_identifier(sim::SimulationParameters)
    return @sprintf("?name=%s?J_pm=%.3f?B=[%.3f,%.3f,%.3f]",
        sim.name,sim.Jpm,sim.B[1],sim.B[2],sim.B[3]) 
end
    
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
#        g1["BZ_grid_density"] = ip.BZ_grid_density

        g = create_group(file, "physical_parameters")
        g["name"] = sim.name
        g["fluxes"] = sim.A
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
    
 
    res =  @showprogress @distributed (vcat) for I = 1:num_K
        k = path.K[I]*0.5
        q = SVector(k[1], k[2], k[3])
        # hard-coded DO g-tensor
        (I, spectral_weight(q, Egrid, sim, ip, g_tensor ) )
    end

    # reorder the data properly
    for (I, x) in res
        Spm[I, :], Spp[I, :], Smagnetic[I, :], bounds[I,:] = x
    end
    
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

    errors=[]

    res = @showprogress @distributed (.+) for _ = 1:ip.n_K_samples  
        q = (1 .- 2 .*(@SVector rand(3)))*4π/8
        p = (1 .- 2 .*(@SVector rand(3)))*4π/8

        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, sim, g_tensor)
		(
        broadened_peaks(S_pm_rs, E_rs, Egrid, ip.broadening_dE ),
		broadened_peaks(S_pp_rs, E_rs, Egrid, ip.broadening_dE ),
        broadened_peaks(S_magnetic_rs::Matrix{Float64}, E_rs, Egrid,
			ip.broadening_dE )
            )
        catch e
            println("Negative dispersion at q=$(q), p=$(p)")
            (Egrid*0, Egrid*0, Egrid*0)
        end
    end
 
    # save the data
    return save_SQW(output_dir, 
        Spm=res[1],
        Spp=res[2],
        Smagnetic=res[3],
        bounds=nothing,
        BZ_path=nothing,
        Egrid=collect(Egrid), sim=sim, ip=ip,
        prefix="integrated"
        )
end


function calc_spinons_along_path(output_dir;
    sim::SimulationParameters,
    path::BZPath)
    bands = @sync @showprogress @distributed (vcat) for k in path.K
        try
            spinon_dispersion(k, sim )[1]'
        catch e
            NaN*zeros(Float64, length(sim.lat.tetra_sites))'
        end
    end

    return save_spinons(output_dir;
        bands=bands,
        BZ_path=path,
        sim=sim)
end


#########################
## PLOTTING

#=
"""
Returns the larges energy eigenvalue it can find (checks the high symmetry
points, and some other random ones)
"""
function estimate_upper_bound(sim)
    test_path = [
    @SVector Float64[0,0,0],
    @SVector Float64[4π/8,0,0],
    @SVector Float64[4π/8,4π/8,0],
    @SVector Float64[2π/8,2π/8,0],
    @SVector Float64[2π/8,2π/8,2π/8]
    ]
    return maximum(
        map(
            k->maximum(spinon_dispersion(k,sim)[1]),
            test_path
        )
    )
end
=#

function plot_spinons(data::Dict)
    d = data["spinon_dispersion"]

    simd = data["physical_parameters"]

    B = simd["B"]

    plot(d["tau"],d["bands"],legend=false,color=:black,lw=0.5)
    xticks!(d["ticks_tau"],d["ticks_label"])
    ylims!(0.,maximum(d["bands"]) )

    bstr = @sprintf("[%.3f,%.3f,%.3f]",B[1],B[2],B[3])  
    if norm( abs.(B/norm(B))- [1,1,1]/√3) < 1e-8
        bstr = @sprintf("%.3f [1,1,1]/\\sqrt{3}", norm(B) )
    elseif norm( abs.(B/norm(B))- [1,1,0]/√2) < 1e-8
        bstr = @sprintf("%.3f [1,1,0]/\\sqrt{2}", norm(B))
    end
    title!(@sprintf("\$J_\\pm=%.3fJ_{yy}, B=%s J_{yy}\$",simd["Jpm"],bstr)  )

    return plot!()
end




""" 
Plots the spectral weight (either the S+S- correlation, or something else)
stored in the simulation stored at the specified path.
"""
function plot_spectral_weight(data::Dict, which="Spm")
    d = data["intensity"]

    p = heatmap(d["tau"],d["W"],real.(d[which])')
    plot!(d["tau"], d["bounds"], linecolor=:white)
    xticks!(d["ticks_tau"],d["ticks_label"])
    plot!(legend=nothing)

    sim =data["physical_parameters"]

    B = sim["B"]
    vline!(d["ticks_tau"], linecolor=:white, linestyle=:dash)
    
    bstr = @sprintf("[%.3f,%.3f,%.3f]",B[1],B[2],B[3])  
    if norm( abs.(B/norm(B))- [1,1,1]/√3) < 1e-8
        bstr = @sprintf("%.3f [1,1,1]/\\sqrt{3}", norm(B) )
    elseif norm( abs.(B/norm(B))- [1,1,0]/√2) < 1e-8
        bstr = @sprintf("%.3f [1,1,0]/\\sqrt{2}", norm(B))
    end
    title!(@sprintf("\$J_\\pm=%.3fJ_{yy}, B=%s J_{yy}\$",sim["Jpm"],bstr))
    return p
end






# some useful defaults

const integration_settings = Dict(
    "very_fast" => IntegrationParameters(n_K_samples=10,   broadening_dE=0.1),
    "fast" =>   IntegrationParameters(n_K_samples=100,  broadening_dE=0.05),
    "slow" =>      IntegrationParameters(n_K_samples=1000, broadening_dE=0.02),
    "very_slow" => IntegrationParameters(n_K_samples=10000,broadening_dE=0.02)
)

