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
