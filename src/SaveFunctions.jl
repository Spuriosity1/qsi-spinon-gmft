using LinearAlgebra
using StaticArrays
using JLD
using HDF5
using Printf
using Optim

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
        g1["version"] = 1.3
#        g1["BZ_grid_density"] = ip.BZ_grid_density

        g = create_group(file, "physical_parameters")
        g["name"] = sim.name
        g["emergent_fluxes"] = calc_fluxes(sim.lat, sim.A)
        g["gauge"] = sim.A
        g["Jpm"] = sim.Jpm
        g["B"] = Vector(sim.B)
        g["L"]=sim.lat.L

        d = create_group(file, "intensity") 
        d["Spm"] = Spm
        d["Spp"] = Spp
        d["Smagnetic"] = Smagnetic
        if bounds !== nothing
            d["bounds"] = bounds
        end
        # a list of K points, such that the I'th S slics corresponds to the I'th K point
        if BZ_path !== nothing
            d["Q_list"] = reduce(hcat, BZ_path.K) 
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
        g1["version"] = 1.3

        g = create_group(file, "physical_parameters")
        g["name"] = sim.name
        g["emergent_fluxes"] = calc_fluxes(sim.lat, sim.A)
        g["gauge"] = sim.A
        g["Jpm"] = sim.Jpm
        g["B"] = Vector(sim.B)
        g["L"]= Vector(sim.lat.L)

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
        g["B"] = Vector(sim.B)
        g["L"]=Vector(sim.lat.L)

        d = create_group(file, "spinon_dispersion") 
        d["bands"] = bands
        # a list of K points, such that the I'th S slics corresponds to the I'th K point
        d["Q_list"] = reduce(hcat, BZ_path.K)
        d["tau"] = BZ_path.t
        d["ticks_tau"] = BZ_path.ticks_t
        d["ticks_label"] = BZ_path.ticks_label
    end
    return name
end


