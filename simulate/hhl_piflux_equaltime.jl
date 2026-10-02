cd(@__DIR__)                       # resolve the ../ paths from the script's dir
using Pkg; Pkg.activate("..")
include("../src/SimFunctions.jl")

using StaticArrays

# ----------------------------------------------------------------------------
# Equal-time (energy-integrated) structure factor in the (hhl) plane for the
# pi-flux configuration (Jpm = 0.2, zero field, 221 cell).
# ----------------------------------------------------------------------------
lat = geom.PyroPrimitive(2, 2, 1)
const sim = SimulationParameters("piflux-hhl",
    lattice = lat,
    A = construct_landau_gauge(lat, [0 0 0 π; 0 π 0 π; 0 0 0 0]),
    Jpm     = 0.2,
    B       = [0.0, 0.0, 0.0],
    n_samples = 10000,
)

# double check phase is right
phi = calc_fluxes(sim)
@assert all(abs.(cos.(phi) .+ 1) .< 1e-6) "A does not correspond to pi-flux"

const csim = CompiledModel(sim)
println("Compiled pi-flux model, spinon mass λ = $(csim.lambda)")

# g-tensor mapping local spin -> magnetic moment; required for Smag (the g-tensor
# / transverse-projector block in corr_at is skipped when g_tensor is nothing).
const G = @SMatrix [0. 0. 0.;
                    0. 0. 0.;
                    1. 0. 0.]

outfile = calc_hhl_equaltime("../output"; csim=csim, N_p=4000, nk=11, hmax=4.0, g_tensor=G)
println("Wrote $(outfile)")
