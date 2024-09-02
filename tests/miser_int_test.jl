struct integral_result
	avg
	var
end

# adapted from https://maxwellrules.com/math/montecarlo_integration.html

@kwdef struct MC_params
	N_points::Int
end

@kwdef struct miser_params 
	MNBS::Int   # Number of points to use for variance exploration
	MNPTS::Int  # Min # points for variance exploration
	pfac::Float64 = 0.1
	alpha::Float64 # splitting fudge factor
end



function MC_integrate(f, xl, xu, N_points::Int64)
	res = MC_kernel(f, xl, xu, N_points::Int64)
	V = prod(xu-xl)
	return integral_result(res.avg*V, res.var*V^2/N_points)
end

function MC_kernel(f, xl, xu, N_points::Int64)
	ndims = length(xl)
	X= rand(ndims, N_points).*(xu - xl) .+ xl
	f_evaluations=f.(eachcol(X))
	μ = sum(f_evaluations)/N_points
	return integral_result(μ, (sum(f_evaluations.*f_evaluations)/N_points - μ^2))
end

function stdev(Y)
	L = length(Y)
	μ = sum(Y)/L
	return sqrt( sum(Y.*Y)/L - μ^2 )
end

function _split_domain(f, xl, xu, N_points::Int)
	# To fix division by zero when variance is zero.
	epsilon = 1e-25

	ndims = length(xu)

	X = rand(ndims, N_points).*(xu - xl) .+ xl
	f_evaluations = f.(eachcol(X))

	# Find midpoints for the bisections
    xmid = 0.5 .* (xl .+ xu)

    # Allocate memory to store variances.
    # axis-1 determines region: 0 -> lower, 1 -> upper
    std_subregions = zeros(Float64, (ndims, 2))
    for dim = 1:ndims
        msk_lu = X[dim, :] .< xmid[dim]

        std_subregions[dim, 1] = stdev(f_evaluations[msk_lu])
        std_subregions[dim, 2] = stdev(f_evaluations[.!msk_lu])
    end

    bisection_dim = argmin(std_subregions[:, 1] + std_subregions[:, 2])

    # If the variance is zero we run into division by zero problems.
    # We plug the hole adding a small epsilon to the variance. The
    # paper proposes a more robust estimator of the variance which is
    # what all "real" implementations use.
    var_l = std_subregions[bisection_dim, 1] .^ 2 .+ epsilon
    var_u = std_subregions[bisection_dim, 2] .^ 2 .+ epsilon

    # Create limits for the split domains that go into
    # the recursion step
    newxu = copy(xu)
    newxl = copy(xl)
    newxu[bisection_dim] = xmid[bisection_dim]
    newxl[bisection_dim] = xmid[bisection_dim]

    return newxu, newxl, var_l, var_u

end

"""
given var_u and var_l, the variances of the upper and lower boxes, distributes
N_points_remaining to minimise variance.
"""
function _allocate_points(var_l, var_u, N_points_remaining, p::miser_params)
    exponent = 1.0 / (p.alpha + 1)

    lower_fraction = var_l^ exponent / (
        var_l^ exponent + var_u^ exponent
    )
    N_points_l = max(ceil(Int64, 
        p.MNPTS + (N_points_remaining - 2 * p.MNPTS) * lower_fraction
		), 1)
	# may borrow a few too many here, this is mild
	N_points_u = max(N_points_remaining - N_points_l, 1)

    return (N_points_l, N_points_u)
end

"""
Integrates f on the domain [xl, xu]
"""
function miser_kernel(f, xl::Vector{Float64}, xu::Vector{Float64}, N_points::Int64, p::miser_params)
	if N_points < p.MNPTS
        return MC_kernel(f, xl, xu, N_points)
    else # recursion path
        # 1) obtain a variance estimate
        npre = ceil(Int64, max(p.pfac * N_points, p.MNPTS))
        # 2) Find best split
        newxu, newxl, var_l, var_u = _split_domain(f, xl, xu, npre)
        # 3) Compute new allocation of points for split regions
        N_points_remaining = N_points - npre
        N_points_l, N_points_u = _allocate_points(
            var_l, var_u, N_points_remaining, p
        )
        # 4) Recurse
        res_l = miser_kernel(f, xl, newxu, N_points_l, p)
        res_u = miser_kernel(f, newxl, xu, N_points_u, p)

        # 5) Combine results

        return integral_result(
            0.5 * (res_l.avg + res_u.avg),
            0.25 * (res_l.var + res_u.var)
        )
    end
end


function miser_integrate(f, xl, xu, N_points::Int64, p::miser_params)
	res = miser_kernel(f, xl, xu, N_points::Int64, p::miser_params)
	V = prod(xu-xl)
	return integral_result(res.avg*V, res.var*V^2/N_points)
end
