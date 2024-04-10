


function lerp(x0, x1, step)
    d = norm(x0-x1)
    T = reverse(range(start=1,stop=0, step=-step/d))
    return T*d, map( t-> t*(x1-x0) + x0, T)
end

function generate_path(point_dict, pathspec, points_per_unit=20)
    #=
    point_dict: A dictionary of the form 'label': [a,b,c]
    pathspec: A list of point labels, e.g. ["\\Gamma", "X", "M", ..., "X"]
    points_per_unit: how densely to make t vectors.

    Returns: tuple,
      t: a 'time' variable
      K: a N*(dim) vector of the lerp'd K points
      t_ticks: the time coordinates for the path points
    =#
    if length(pathspec) == 0
        print("Path has length zero, exiting")
        return [], [], []
    end

    
    
    t_ticks = [0.0]
    T = [0.0]
    p0 = point_dict[pathspec[1]]
    K = [p0]
    step = 1/points_per_unit
    
    for label in pathspec[2:length(pathspec)]
        p = point_dict[label]
        tau, K_loc = lerp(p0, p, step)
        append!(T, collect(tau).+last(T))
        push!(t_ticks, last(T) )
        append!(K, K_loc)
        
        p0 = p
    end

    return T, K, t_ticks
end
