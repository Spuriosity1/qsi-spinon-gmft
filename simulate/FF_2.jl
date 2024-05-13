using Distributed

include("driver.jl")

A_FF =  load_A("gaugefiles/FF_even_pi18.gauge")

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

Returns the (111) field needed to realise the prescribed field
"""
function FF_111_B(desired_Φ0, Jpm)
    # x = Φ0 -> ( 1 - 2*cos((π/2 - Φ0)*2/3) )^-1;
    xp = g0_g123_ratio_from_flux(desired_Φ0)
    
    B = sqrt(Jpm* 6*(1-xp)/ (xp*5/9 - 5) )
    return B
end


Φ0 = π/3

simlist = map(
    Jpm->SimulationParameters("FF",
        A=A_FF, 
        Jpm=Jpm,
        B=FF_111_B(Φ0, Jpm)* [1.,1.,1.]/sqrt(3),
        nsample=1000,
        kappa=2.0
    ),
    [-0.01, -0.02, -0.05, -0.1]
)

for (i, sim) in enumerate(simlist)
    @printf("Running simulation %d of %d\n", i, length(simlist))
    run_sim(
        data_dir="output/",
        figure_dir="figures/",
        sim=sim, 
        integral_params=integration_settings["very_fast"],
        k_density_spinon_dispersion=20,
        k_density_specweight=10
        )
end


    
