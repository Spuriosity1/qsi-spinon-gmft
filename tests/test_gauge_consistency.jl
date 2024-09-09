include("../src/SpinonStructure.jl")
include("../src/BZMath.jl")


using .BZmath
using .SpinonStructure
using StaticArrays
using SparseArrays
using LinearAlgebra
using BenchmarkTools
using ProgressMeter

# number of test points to check
const NCHECK=100

# const integration_settings = Dict(
#     "very_fast" => IntegrationParameters(n_K_samples=10,   broadening_dE=0.1),
#     "fast" =>   IntegrationParameters(n_K_samples=100,  broadening_dE=0.05),
#     "slow" =>      IntegrationParameters(n_K_samples=1000, broadening_dE=0.02),
#     "very_slow" => IntegrationParameters(n_K_samples=10000,broadening_dE=0.02),
#     "ultra_slow" => IntegrationParameters(n_K_samples=100000,broadening_dE=0.02)
# )

lat = geom.PyroPrimitive(2,3,4)

bfield = [0.1,0.4,0.9]
#bfield = [0.,0.,0.]

# the control / reference
sim1 = SimulationParameters("0flux_std_gauge",
    lattice=lat,
    A=zeros(div(length(lat.tetra_sites),2),4),
    Jpm=-0.046,
    B=bfield,
    n_samples=10000
    )

ch1 = SpinonStructure.CompiledHamiltonian(sim1)

function generate_gauged_sim(Γ, shiftg)

    sim_test = SimulationParameters("0flux_random",
    lattice=lat,
    # A=[1 1 0 0; 0 0 -1 -1; 1 1 0 0; 1 1 0 0],
    A=lattice_gradient(lat, Γ) .+ shiftg,
    Jpm=-0.046,
    B=bfield,
    n_samples=10000
    )

    @assert norm(calc_fluxes(sim1)-calc_fluxes(sim_test)) < 1e-10

    return SpinonStructure.CompiledHamiltonian(sim_test)
end


n_failures = 0;



println("Testing M's gauge convariance for local gauges... (0-flux)")
@showprogress for n=1:NCHECK  
    Q = SVector{3}(rand(3))*2π;
    H1 = SpinonStructure.calc_hopping(ch1, Q)

    Γ = 2π .*rand(Float64, size(lat.tetra_sites)).-π


    ch2 = generate_gauged_sim(Γ, [0 0 0 0])
    
    H2 = SpinonStructure.calc_hopping(ch2, Q)

    gauge = diagm(exp.(+1im.*Γ))
    diff = gauge'*H2*gauge  - H1
    if norm(diff) > 1e-10
        println("TEST FAILED: gauge transform not gauge")
        display(norm.(diff))
        global n_failures += 1;
    end
end



println("Testing M's gauge covariance for random large gauges...")
@showprogress for n=1:NCHECK  
    Q = SVector{3}(rand(3))*2π;

    # the sublat gauge (make divergence free - divergenceful changes
	# correspond to a different local gauge
	σ = (rand(4).-0.5)*4π
	lam = sum(σ)/4
	σ .-= lam

	Bmat = reduce(vcat,[2*eμ' for eμ in geom.pyro])
		deltaQ = Bmat[1:3,:] \ (σ[1:3])


    gaugevec = zeros( size(lat.tetra_sites))
    ch2 = generate_gauged_sim(gaugevec, σ')
 
#    H1 = SpinonStructure.calc_hopping(ch1, Q)
    #H2 = SpinonStructure.calc_hopping(ch2, Q+deltaQ)
    

	H2 = SpinonStructure.calc_hopping(ch2, Q)
    H1 = SpinonStructure.calc_hopping(ch1, Q-deltaQ)


	gauge = diagm(exp.(+1im.*gaugevec))
	diff =  H2  - gauge*H1*gauge'
	if norm(diff) > 1e-10
        println("TEST FAILED: gauge transform not gauge")
        display(norm.(diff))
        global n_failures += 1
    end
end






if n_failures > 0
    println("$(n_failures)/$(NCHECK) failed.");
    exit(n_failures)
end

