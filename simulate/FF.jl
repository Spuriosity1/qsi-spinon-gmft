using Distributed

include("driver-headless.jl")


const al1 = [ 0 0 -π/2 π/2; 0 0 0 0; 0 0 π/2 0 ]


"""
Returns the ratio of g0 to g1 in the special case B || [111], given a particular value of 
Φ0, the flux on plaquette 0 
The other three plaquettes have flux -Φ0/3
"""
g0_g123_ratio_from_flux(desired_Φ0) = (1-4*cos(desired_Φ0/3)^2)^-1

"""
FF_111_B(desired_Φ, Jpm)
-> B

x = g0/g1 <0

Returns the (111) field needed to realise the prescribed flux
"""
function FF_111_B(desired_Φ0, Jpm)
    # x = Φ0 -> ( 1 - 2*cos((π/2 - Φ0)*2/3) )^-1;
    xp = g0_g123_ratio_from_flux(desired_Φ0)
    
    B = sqrt(Jpm* 6*(1-xp)/ (xp*5/9 - 5) )
    return B
end




function sim_factory(N::Int, Jpm::Float64)
    lat = geom.PyroPrimitive(4*N,1,4*N)

    A_FF = construct_landau_gauge(lat, al1/N)

    flux_state=calc_fluxes(lat, A_FF)

    mean_fluxes = sum.(eachcol(flux_state))/size(flux_state,1)

    # check homogeneity
    max_err = sqrt(maximum((flux_state .- mean_fluxes').^2))
    println("Largest deviation from mean flux is ", max_err)


    Φ0 = abs(mean_fluxes[1])
    println("$(N) -> phi = $(Φ0)")

    return SimulationParameters("FF",
        A=A_FF, 
        lattice=lat,
        Jpm=Jpm,
        B=FF_111_B(Φ0, Jpm)* [1.,1.,1.]/sqrt(3),
        n_samples=1000,
        kappa=2.0
        )
end


Jpm = -0.04

simlist = map(
N->sim_factory(N, Jpm),
    [2,3,4]
)


for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output/",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["fast"],
        k_density_spinon_dispersion=20,
        k_density_specweight=5
        )
end


    
