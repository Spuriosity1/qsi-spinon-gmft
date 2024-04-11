

function lerp(x0, x1, step)
    d = norm(x0-x1)
    T = reverse(range(start=1,stop=0, step=-step/d))
    return T*d, map( t-> t*(x1-x0) + x0, T)
end



"""
essentially a named tuple
"""
struct BZPath
    t::Vector{Float64}
    K::Vector{Vector{Float64}}
    ticks_t::Vector{Float64}
    ticks_label::Vector{String}
end

"""
    generate_path(point_dict, pathspec, points_per_unit=20)

Arguments:
    point_dict: A dictionary of the form 'label': [a,b,c]
    pathspec: A list of point labels, e.g. ["\\Gamma", "X", "M", ..., "X"]
    points_per_unit: how densely to make t vectors.

Returns:    
    a BZpath object
    
"""
function generate_path(point_dict, pathspec; points_per_unit=20, K_units=1.)

    if length(pathspec) == 0
        print("Path has length zero, exiting")
        return [], [], []
    end

    
    
    t_ticks = [0.0]
    T = [0.0]
    p0 = point_dict[pathspec[1]]*K_units
    K = [p0]
    step = 1/points_per_unit
    
    for label in pathspec[2:length(pathspec)]
        p = point_dict[label]*K_units
        tau, K_loc = lerp(p0, p, step*K_units)
        append!(T, collect(tau).+last(T))
        push!(t_ticks, last(T) )
        append!(K, K_loc)
        
        p0 = p
    end
    
    return BZPath(T, K, t_ticks, pathspec)
end


