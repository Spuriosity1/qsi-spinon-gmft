include("SpinonStructure.jl")
include("BZMath.jl")
using .BZmath
using .SpinonStructure
using LinearAlgebra
using StaticArrays
using ProgressMeter
using Plots
using JLD
using HDF5
using Printf


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
        S::Matrix{ComplexF64},
        bounds::Matrix{Float64},
        Egrid::Vector{Float64},
        BZ_path::BZPath,
        sim::SimulationParameters, 
        ip::IntegrationParameters,
        prefix="SQW")

    name = output_dir*"/"*prefix*sim_identifier(sim)
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
        d["S"] = S
        d["bounds"] = bounds
        # a list of K points, such that the I'th S slics corresponds to the I'th K point
        d["Q_list"] = BZ_path.K 
        d["tau"] = BZ_path.t
        d["ticks_tau"] = BZ_path.ticks_t
        d["ticks_label"] = BZ_path.ticks_label
        d["W"] = Egrid
    end
    return name
end





"""
Plots integrated spectral weight over the whole Brillouin zone
"""
function integrated_specweight(sim::SimulationParameters, 
						 integral_params::IntegrationParameters,
        Egrid::Vector{Float64}
						 )
    Sω = zeros(ComplexF64,size(Egrid))

    p = Progress(integral_params.n_K_samples)
    @Threads.threads for _ = 1:integral_params.n_K_samples  
        q = (1 .- 2 .*(@SVector rand(3)))*4π/8
        k = (1 .- 2 .*(@SVector rand(3)))*4π/8
        E, S = SpinonStructure.specweight_at(q, k, sim)


        Sω += map( 
            e-> sum( 
                [ S*Lorentzian(e - E, integral_params.broadening_dE)
                for (E,S) in zip(Enm,Snm)
                    ]), 
            Egrid
            )
        next!(p)
    end
    finish!(p)
    return Sω
end  


###########
# expensive calculations
"""
calc_spectral_weight(sim::SimulationParameters,
    ip::IntegrationParameters, Egrid::Vector{Float64}, path::BZPath)

    @param sim the simulation parameters
    @param ip the numerical constants, passed directly to spectral_weight
    @param Egrid the grid for energy values
    @param path the path through the Brillouin zone

    Saves the results to
"""
function calc_spectral_weight_along_path(sim::SimulationParameters,
    ip::IntegrationParameters, Egrid::Vector{Float64}, path::BZPath,
    output_dir::String)
    
    num_K = length(path.K)
    
    S = zeros(ComplexF64, num_K, length(Egrid))
    bounds = zeros(Float64, num_K, 2)
    
    p = Progress(num_K)
 
    Threads.@threads for I = 1:num_K
        k = path.K[I]*0.5
        q = SVector(k[1], k[2], k[3])
        S[I, :], bounds[I,:] = spectral_weight(q, Egrid, sim, ip)
        next!(p)
    end
    
    finish!(p)

    # save the data
    return save_SQW(output_dir, S=S, bounds=bounds, BZ_path=path,
        Egrid=collect(Egrid), sim=sim, ip=ip)
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

function plot_spinons(sim::SimulationParameters, path::BZPath)
    E = []
    @showprogress for k in path.K
        push!(E, spinon_dispersion(k, sim )[1]')
    end
    E = reduce(vcat, E)

    plot(path.t,E,legend=false,color=:black,lw=0.5)
    xticks!(path.ticks_t, path.ticks_label)
    ylims!(0.,maximum(E))

    bstr = @sprintf("[%.3f,%.3f,%.3f]",sim.B[1],sim.B[2],sim.B[3])  
    if norm( abs.(sim.B/norm(sim.B))- [1,1,1]/√3) < 1e-8
        bstr = @sprintf("%.3f [1,1,1]/\\sqrt{3}", norm(sim.B) )
    elseif norm( abs.(sim.B/norm(sim.B))- [1,1,0]/√2) < 1e-8
        bstr = @sprintf("%.3f [1,1,0]/\\sqrt{2}", norm(sim.B))
    end
    title!(@sprintf("\$J_\\pm=%.3fJ_{yy}, B=%s J_{yy}\$",sim.Jpm,bstr)  )

    return plot!()
end


""" 
Plots the spectral weight (either the S+S- correlation, or something else)
stored in the simulation stored at the specified path.
"""
function plot_spectral_weight(data::Dict)
    d = data["intensity"]
    p = heatmap(d["tau"],d["W"],real.(d["S"])')
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
