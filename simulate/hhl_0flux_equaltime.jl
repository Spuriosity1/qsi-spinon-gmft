cd(@__DIR__)                       # resolve the ../ paths from the script's dir
using Pkg; Pkg.activate("..")
include("../src/SimFunctions.jl")

using StaticArrays

# ----------------------------------------------------------------------------
# Equal-time (energy-integrated) structure factor in the (hhl) plane for the
# 0-flux configuration (Jpm = -0.04, zero field, smallest primitive cell).
# ----------------------------------------------------------------------------
const sim = SimulationParameters("0flux-hhl",
    lattice = geom.PyroPrimitive(1, 1, 1),
    A       = [0 0 0 0],
    Jpm     = -0.04,
    B       = [0.0, 0.0, 0.0],
    n_samples = 10000,
)

# sanity check: this really is the 0-flux ground state for these couplings
@assert all(abs.(calc_fluxes(sim)) .< 1e-9) "A does not correspond to 0-flux"

const csim = CompiledModel(sim)
println("Compiled 0-flux model, spinon mass λ = $(csim.lambda)")

outfile = calc_hhl_equaltime("../output"; csim=csim, N_p=40000, nk=121, hmax=4.0)
println("Wrote $(outfile)")
