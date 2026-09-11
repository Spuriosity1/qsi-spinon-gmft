using ProgressMeter
include("SaveFunctions.jl")

using StaticArrays


"""
    equaltime_corr(q, csim, N_p; g_tensor=nothing)

Energy-integrated structure factor at a single external momentum `q`, MC-averaged
over the internal spinon momentum `p` sampled uniformly over one primitive
reciprocal cell (which tiles the BZ).

Each spectral peak of weight `S` sits under `∫dω = S`, so integrating over all ω
returns the sum of the raw `corr_at` weights -- no energy grid or broadening
required. Returns `(Spm, Spp, Smag)`; any component is `NaN` if every sample hit
a negative-dispersion `DomainError`.
"""
function equaltime_corr(q::SVector{3,Float64}, csim::CompiledModel, N_p::Int;
        g_tensor=nothing)
    acc_pm = 0.0
    acc_pp = 0.0
    acc_mag = 0.0
    ngood = 0
    for _ in 1:N_p
        p = geom.primitive_recip_basis * (SVector{3}(rand(), rand(), rand()) .- 0.5)
        try
            _, S_pm, S_pp, S_mag = corr_at(q, p, csim, g_tensor)
            acc_pm  += real(sum(S_pm))
            acc_pp  += real(sum(S_pp))
            acc_mag += real(sum(S_mag))
            ngood += 1
        catch e
            e isa DomainError || rethrow(e)   # skip negative-dispersion points
        end
    end
    ngood == 0 && return (NaN, NaN, NaN)
    return (acc_pm / ngood, acc_pp / ngood, acc_mag / ngood)
end


"""
    calc_hhl_equaltime(output_dir; csim, N_p=4000, nk=121, hmax=4.0, g_tensor=nothing)

Sweep the (hhl) plane, evaluating the equal-time structure factor at each
`q = (2π/8)·(h, h, l)` with `h, l ∈ [-hmax, hmax]` r.l.u. on an `nk × nk` grid.
`N_p` is the number of internal-momentum MC samples per q. Saves Spm/Spp/Smag to
an HDF5 file and returns its path.

Reciprocal-space units: the integer lattice uses a conventional cubic cell of
side 8, so a reciprocal-lattice unit (r.l.u.) is 2π/8.
"""
function calc_hhl_equaltime(output_dir::String;
        csim::CompiledModel,
        N_p::Int=4000,
        nk::Int=121,
        hmax::Float64=4.0,
        g_tensor=nothing)

    rlu   = 2π / 8
    haxis = collect(range(-hmax, hmax, length=nk))
    laxis = collect(range(-hmax, hmax, length=nk))

    Spm  = fill(NaN, nk, nk)       # indexed [il, ih]
    Spp  = fill(NaN, nk, nk)
    Smag = fill(NaN, nk, nk)

    prog = Progress(nk * nk, desc="(hhl) equal-time: ")
    Threads.@threads for ih in 1:nk
        h = haxis[ih]
        for il in 1:nk
            l = laxis[il]
            q = SVector{3,Float64}(rlu * h, rlu * h, rlu * l)
            Spm[il, ih], Spp[il, ih], Smag[il, ih] =
                equaltime_corr(q, csim, N_p; g_tensor=g_tensor)
            next!(prog)
        end
    end
    finish!(prog)

    return save_hhl_equaltime(output_dir;
        haxis=haxis, laxis=laxis, Spm=Spm, Spp=Spp, Smag=Smag,
        csim=csim, N_p=N_p)
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
    sim::CompiledModel,
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
        k = path.K[I]
        q = SVector(k[1], k[2], k[3])

        try
            tmp_intensity = spectral_weight(q, Egrid, sim, ip, g_tensor )
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
        sim=sim.sim, 
        ip=ip
        )
end



function calc_integrated_specweight(
    output_dir::String;
    csim::CompiledModel,
    ip::IntegrationParameters, Egrid::Vector{Float64},g_tensor=nothing)


#=
    Sω_pm = zeros(ComplexF64,size(Egrid))
    Sω_pp = zeros(ComplexF64,size(Egrid))
    Sω_magnetic = zeros(Float64,size(Egrid))
=#

    prog = Progress(ip.n_K_samples, desc="BZ Average")

    samples_per_trial = ip.n_K_samples ÷ Threads.nthreads()

    deltaE = maximum(Egrid[2:end] - Egrid[1:end-1])

    
    function single_trial()        
        Spm_res = zeros(ComplexF64, length(Egrid) ) 
        Spp_res = zeros(ComplexF64, length(Egrid) ) 
        Smagnetic_res = zeros(Float64, length(Egrid) ) 


        for _ in 1:samples_per_trial
            q = (1 .- 2 .*(@SVector rand(3)))*4π/8
            p = (1 .- 2 .*(@SVector rand(3)))*4π/8

            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, csim, g_tensor)

            broadened_peaks!(Spm_res, S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, deltaE )
            broadened_peaks!(Spp_res, S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, deltaE )
            broadened_peaks!(Smagnetic_res, S_magnetic_rs::Matrix{Float64}, E_rs, Egrid, deltaE )

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
                sim=csim.sim, ip=ip, prefix="integrated")
end




# single-threaded version
function calc_integrated_S(;
    csim::CompiledModel,
    ip::IntegrationParameters,
    Egrid::Vector{Float64},
    g_tensor=nothing
    )
 
    Spm_res = zeros(ComplexF64, length(Egrid) ) 
    Spp_res = zeros(ComplexF64, length(Egrid) ) 
    Smagnetic_res = zeros(Float64, length(Egrid) ) 


    deltaE = maximum(Egrid[2:end] - Egrid[1:end-1])

    for _ = 1:ip.n_K_samples  
        q = (1 .- 2 .*(@SVector rand(3)))*4π/8
        p = (1 .- 2 .*(@SVector rand(3)))*4π/8

        try
            E_rs, S_pm_rs, S_pp_rs, S_magnetic_rs = corr_at(q, p, csim, g_tensor) 
            Spm_res += broadened_peaks(S_pm_rs::Matrix{ComplexF64}, E_rs, Egrid, deltaE )
            Spp_res += broadened_peaks(S_pp_rs::Matrix{ComplexF64}, E_rs, Egrid, deltaE )
            Smagnetic_res += broadened_peaks(S_magnetic_rs::Matrix{Float64}, E_rs, Egrid,
                                         deltaE )
                    
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
	Spm_res, Spp_res, Smagnetic_res = calc_integrated_S(                                                          csim=CompiledModel(this_sim),
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
    csim::CompiledModel,
    path::BZPath,
    kshift=[0.,0.,0.]
    )
    num_K = length(path.K)
    bands = zeros(Float64, num_K, length(csim.sim.lat.tetra_sites))
    
    prog=Progress(num_K, desc="Spinon Dispersion: ") 
    @Threads.threads for J=1:num_K
        bands[J,:] = begin
            try
                spinon_dispersion(path.K[J]+kshift, csim)[1]'
            catch e
                if e isa DomainError
                    NaN*ones(Float64, 1,length(csim.sim.lat.tetra_sites))
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
        sim=csim.sim)
end

"""
Jpm -> The Jpm parameter
B -> a (3,1) vector of magnetic int_fields
Returns a 4-component Vector{Float64} of g-values,
corresponding to plaquettes with normals in the directions set by geom.pyro
i.e. [1,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1] respectively
"""
function Jring(Jpm::Float64, B::Union{Vector{Float64},SVector{3,Float64}})
    # assumes units of Jy
    ring_normals =  [1. 1 1; 1 -1 -1;-1 1 -1;-1 -1 1]

    return (3 * Jpm^3/2) .- (5/12) .* Jpm^2 .* (ring_normals * B).^2
end

"""
Calculates the minimum flux of the plaquettes for the hamiltonian
        H = ∑ g_μ cos( Φ_μ )
returns [Φ1,Φ2,Φ3,Φ4]
"""
function optimal_flux(g)
    objective(x) = g[1]*cos(-sum(x)) + sum(g[2:4].*cos.(x))
    function gradient!(G,x)
        G[1] = g[1]*sin(sum(x)) - g[2]*sin(x[1])
        G[2] = g[1]*sin(sum(x)) - g[3]*sin(x[2])
        G[3] = g[1]*sin(sum(x)) - g[4]*sin(x[3])
    end
    lower = [0.,0.,0.].-1e-5
    upper = [π,π,π]

    results = [
        optimize(objective, lower, upper, x0)
        for x0 in [zeros(3), (π-0.1)*ones(3), π/3*ones(3), 2π/3*ones(3) ]
    ]
    I = argmin(res.minimum for res in results)
    min_x = results[I].minimizer
    X= mod.([-sum(min_x), min_x...].+π,2π).-π
    if X[1] > 0
        X .*= -1
    end
    return X
end

        # some useful defaults

const integration_settings = Dict(
    "very_fast" =>  IntegrationParameters(n_K_samples=10    ),
    "fast" =>       IntegrationParameters(n_K_samples=100   ),
    "slow" =>       IntegrationParameters(n_K_samples=1000  ),
    "very_slow" =>  IntegrationParameters(n_K_samples=10000 ),
    "ultra_slow" => IntegrationParameters(n_K_samples=100000)
)

