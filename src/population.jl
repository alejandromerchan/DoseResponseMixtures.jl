"""
    AbstractPopulation

Supertype for all population models in DoseResponseMixtures. Each concrete
subtype represents a different model of how individual tolerance is
distributed in a population.
"""
abstract type AbstractPopulation end

"""
    UnimodalPopulation(mean_log_lc50, sd_log_lc50, slope)

A population whose log-LC50 is normally distributed with mean
`mean_log_lc50` and standard deviation `sd_log_lc50`. The `slope` parameter
controls the steepness of the dose-response curve at the individual level.

# Arguments
- `mean_log_lc50::Float64`: Mean of the log-LC50 distribution
- `sd_log_lc50::Float64`: Standard deviation of log-LC50 (≥ 0; use 0 for
  a homogeneous population with fixed tolerance)
- `slope::Float64`: Slope of the individual-level log-logistic response (> 0)
"""
struct UnimodalPopulation <: AbstractPopulation
    mean_log_lc50::Float64
    sd_log_lc50::Float64
    slope::Float64

    function UnimodalPopulation(mean_log_lc50, sd_log_lc50, slope)
        sd_log_lc50 >= 0 || throw(ArgumentError("sd_log_lc50 must be ≥ 0"))
        slope > 0 || throw(ArgumentError("slope must be > 0"))
        return new(Float64(mean_log_lc50), Float64(sd_log_lc50), Float64(slope))
    end
end

"""
    HWEMixturePopulation(allele_frequency, mean_log_lc50_susceptible,
                         resistance_magnitude, dominance, sd_log_lc50, slope)

A three-component mixture population representing a single biallelic
major-effect resistance locus under Hardy-Weinberg equilibrium.

Under HWE with R allele frequency `q = allele_frequency`, the three
genotype classes have frequencies `(1-q)²` (SS), `2q(1-q)` (SR), and `q²`
(RR). Each genotype has a normal log-LC50 distribution with mean determined
by the dominance coefficient and shared standard deviation `sd_log_lc50`:

  - μ_SS = mean_log_lc50_susceptible
  - μ_SR = mean_log_lc50_susceptible + dominance × resistance_magnitude
  - μ_RR = mean_log_lc50_susceptible + resistance_magnitude

# Arguments
- `allele_frequency::Float64`: R allele frequency q ∈ [0, 1]
- `mean_log_lc50_susceptible::Float64`: Mean log-LC50 of the SS class
- `resistance_magnitude::Float64`: Log-scale tolerance difference between
  RR and SS (any real number; positive means R confers higher tolerance)
- `dominance::Float64`: Dominance coefficient h ∈ [0, 1]; 0 = fully
  recessive R, 1 = fully dominant R, 0.5 = additive
- `sd_log_lc50::Float64`: Shared within-genotype standard deviation (≥ 0)
- `slope::Float64`: Individual-level log-logistic slope (> 0)

# Notes
The case `resistance_magnitude == 0` collapses the model to a single
normal centered at `mean_log_lc50_susceptible`, regardless of `q` or `h`.

Overdominance (h > 1) and underdominance (h < 0) are not supported by
this constructor; if needed they can be added in a separate type.
"""
struct HWEMixturePopulation <: AbstractPopulation
    allele_frequency::Float64
    mean_log_lc50_susceptible::Float64
    resistance_magnitude::Float64
    dominance::Float64
    sd_log_lc50::Float64
    slope::Float64

    function HWEMixturePopulation(
        allele_frequency,
        mean_log_lc50_susceptible,
        resistance_magnitude,
        dominance,
        sd_log_lc50,
        slope,
    )
        0 <= allele_frequency <= 1 ||
            throw(ArgumentError("allele_frequency must be in [0, 1]"))
        0 <= dominance <= 1 || throw(ArgumentError("dominance must be in [0, 1]"))
        sd_log_lc50 >= 0 || throw(ArgumentError("sd_log_lc50 must be ≥ 0"))
        slope > 0 || throw(ArgumentError("slope must be > 0"))
        return new(
            Float64(allele_frequency),
            Float64(mean_log_lc50_susceptible),
            Float64(resistance_magnitude),
            Float64(dominance),
            Float64(sd_log_lc50),
            Float64(slope),
        )
    end
end

"""
    genotype_means(pop::HWEMixturePopulation) -> NamedTuple

Return the three genotype-specific mean log-LC50 values as a NamedTuple
with fields `SS`, `SR`, `RR`.
"""
function genotype_means(pop::HWEMixturePopulation)
    μ_SS = pop.mean_log_lc50_susceptible
    μ_RR = μ_SS + pop.resistance_magnitude
    μ_SR = μ_SS + pop.dominance * pop.resistance_magnitude
    return (SS = μ_SS, SR = μ_SR, RR = μ_RR)
end

"""
    genotype_frequencies(pop::HWEMixturePopulation) -> NamedTuple

Return the Hardy-Weinberg genotype frequencies as a NamedTuple with
fields `SS`, `SR`, `RR`.
"""
function genotype_frequencies(pop::HWEMixturePopulation)
    q = pop.allele_frequency
    return (SS = (1 - q)^2, SR = 2 * q * (1 - q), RR = q^2)
end
