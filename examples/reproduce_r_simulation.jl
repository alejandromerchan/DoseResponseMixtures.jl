# Corresponds to Simulation_Parallel.R in the original R code.
# Validates that the Julia generative model + 2PL estimator reproduces the
# distribution of LC50 estimates from the original R simulation.
#
# Before running:
#   using Pkg
#   Pkg.activate("..")
#   Pkg.instantiate()
#
# StableRNGs must be available (install from the test environment or add it
# to the parent project: Pkg.add("StableRNGs")).

using DoseResponseMixtures
using StableRNGs: StableRNG
using Statistics: mean, std, quantile

# --- Simulation parameters (matching original R script) ---
const TRUE_LC50 = 5.0
const TRUE_SLOPE = 10.0
const SD_LOG_LC50 = 0.2         # approximate equivalent of R's tolerance spread
const DOSES = 1:10
const N_REPLICATES = 5
const N_PER_REPLICATE = 1000
const N_TRIALS = 1000

pop = UnimodalPopulation(log(TRUE_LC50), SD_LOG_LC50, TRUE_SLOPE)
rng = StableRNG(42)

lc50_estimates = Vector{Float64}(undef, N_TRIALS)

for trial = 1:N_TRIALS
    df = simulate_bioassay(pop, DOSES, N_REPLICATES, N_PER_REPLICATE; rng)
    result = fit_2pl(df)
    lc50_estimates[trial] = result.lc50
end

println("LC50 estimate summary over $N_TRIALS trials:")
println("  True LC50 : $TRUE_LC50")
println("  Mean      : $(round(mean(lc50_estimates); digits=3))")
println("  SD        : $(round(std(lc50_estimates); digits=3))")
println("  2.5%      : $(round(quantile(lc50_estimates, 0.025); digits=3))")
println("  50%       : $(round(quantile(lc50_estimates, 0.50); digits=3))")
println("  97.5%     : $(round(quantile(lc50_estimates, 0.975); digits=3))")

# Optionally write results to CSV for plotting in a later milestone
results_path = joinpath(@__DIR__, "lc50_estimates.csv")
open(results_path, "w") do io
    println(io, "lc50")
    for x in lc50_estimates
        println(io, x)
    end
end
println("\nResults written to: $results_path")
