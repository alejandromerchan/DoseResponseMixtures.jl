# Explore HWEMixturePopulation: how the log-tolerance distribution morphs
# with allele frequency q and dominance h.
#
# Before running:
#   using Pkg
#   Pkg.activate("..")
#   Pkg.instantiate()
#
# StableRNGs and Plots must be available (install from the test environment
# or add to the parent project).

using Revise
using DoseResponseMixtures
using StableRNGs: StableRNG
using Plots

const N_SAMPLES = 50_000
const BASE_MU = log(5.0)
const RESIST_MAG = 2.0
const SD = 0.4
const SLOPE = 10.0

# --- Figure 1: vary q at fixed h = 0.5 ---

q_values = [0.05, 0.25, 0.5, 0.75]
h_fixed = 0.5

fig1 = plot(layout = (2, 2), size = (900, 700), title_location = :left)

for (i, q) in enumerate(q_values)
    pop = HWEMixturePopulation(q, BASE_MU, RESIST_MAG, h_fixed, SD, SLOPE)
    rng = StableRNG(100 + i)
    samples = sample_log_tolerances(pop, N_SAMPLES; rng)
    means = genotype_means(pop)
    freqs = genotype_frequencies(pop)

    println("q = $q, h = $h_fixed")
    println(
        "  Frequencies: SS=$(round(freqs.SS; digits=4))  SR=$(round(freqs.SR; digits=4))  RR=$(round(freqs.RR; digits=4))",
    )
    println(
        "  Means:       SS=$(round(means.SS; digits=3))  SR=$(round(means.SR; digits=3))  RR=$(round(means.RR; digits=3))",
    )

    histogram!(
        fig1,
        samples;
        subplot = i,
        bins = 80,
        normalize = :pdf,
        label = "",
        title = "q = $q",
        xlabel = "log-LC50",
        ylabel = "density",
        alpha = 0.7,
        color = :steelblue,
    )
    vline!(fig1, [means.SS]; subplot = i, label = "μ_SS", color = :green, lw = 2)
    vline!(fig1, [means.SR]; subplot = i, label = "μ_SR", color = :orange, lw = 2)
    vline!(fig1, [means.RR]; subplot = i, label = "μ_RR", color = :red, lw = 2)
end

plot!(fig1; plot_title = "HWE mixture: varying q (h = $h_fixed)")
savefig(fig1, joinpath(@__DIR__, "hwe_vary_q.png"))
println("\nFigure 1 saved to examples/hwe_vary_q.png")

# --- Figure 2: vary h at fixed q = 0.5 ---

h_values = [0.0, 0.33, 0.66, 1.0]
q_fixed = 0.5

fig2 = plot(layout = (2, 2), size = (900, 700), title_location = :left)

println()
for (i, h) in enumerate(h_values)
    pop = HWEMixturePopulation(q_fixed, BASE_MU, RESIST_MAG, h, SD, SLOPE)
    rng = StableRNG(200 + i)
    samples = sample_log_tolerances(pop, N_SAMPLES; rng)
    means = genotype_means(pop)
    freqs = genotype_frequencies(pop)

    println("q = $q_fixed, h = $h")
    println(
        "  Frequencies: SS=$(round(freqs.SS; digits=4))  SR=$(round(freqs.SR; digits=4))  RR=$(round(freqs.RR; digits=4))",
    )
    println(
        "  Means:       SS=$(round(means.SS; digits=3))  SR=$(round(means.SR; digits=3))  RR=$(round(means.RR; digits=3))",
    )

    histogram!(
        fig2,
        samples;
        subplot = i,
        bins = 80,
        normalize = :pdf,
        label = "",
        title = "h = $h",
        xlabel = "log-LC50",
        ylabel = "density",
        alpha = 0.7,
        color = :orchid,
    )
    vline!(fig2, [means.SS]; subplot = i, label = "μ_SS", color = :green, lw = 2)
    vline!(fig2, [means.SR]; subplot = i, label = "μ_SR", color = :orange, lw = 2)
    vline!(fig2, [means.RR]; subplot = i, label = "μ_RR", color = :red, lw = 2)
end

plot!(fig2; plot_title = "HWE mixture: varying h (q = $q_fixed)")
savefig(fig2, joinpath(@__DIR__, "hwe_vary_h.png"))
println("\nFigure 2 saved to examples/hwe_vary_h.png")
