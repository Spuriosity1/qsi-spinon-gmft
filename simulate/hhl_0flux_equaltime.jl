#!/usr/bin/env -S julia --threads=auto --project=..
#SBATCH --job-name=hhl-0flux
#SBATCH --output=hhl-0flux-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --mem=8G
#
# Submit directly:  sbatch hhl_0flux_equaltime.jl
# `--threads=auto` lets Julia use the allocated CPUs; on a shared node pin it
# instead with `export JULIA_NUM_THREADS=$SLURM_CPUS_PER_TASK`.

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

outfile = calc_hhl_equaltime("../output"; csim=csim, N_p=4000, nk=121, hmax=4.0)
println("Wrote $(outfile)")
