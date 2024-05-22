include("src/SpinonStructure.jl")
using .SpinonStructure
using JLD
using HDF5
using StaticArrays

files = readdir("output")

for f in files
    println(f)
    h5open("output/"*f, "r+") do file
        println(typeof(file))
        A = read(file, "physical_parameters/fluxes")

        file["physical_parameters"]["emergent_fluxes"] = calc_fluxes( A )
        file["physical_parameters"]["gauge"] =  A
        delete_object(file["physical_parameters"], "fluxes")
    end
#    save("output_fixed/"*f, file)
end

