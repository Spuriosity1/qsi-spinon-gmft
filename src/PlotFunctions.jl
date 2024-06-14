
using Plots
using JLD
using HDF5
using Printf


include("BZMath.jl")
#include("SpinonStructure.jl")

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




