using DataFrames: DataFrame
using GLM: glm, coef, Binomial, LogitLink, FrequencyWeights
using StatsModels: @formula

"""
    fit_2pl(data::DataFrame) -> NamedTuple

Fit a classical two-parameter log-logistic (2PL) dose-response model to
bioassay data using a binomial GLM with logit link.

Assumes a single unimodal tolerance distribution. The model is:

    logit(p) = β₀ + β₁ · log_dose

where `p` is the probability of mortality at a given log-dose. LC50 is
recovered from the GLM coefficients as:

    LC50 = exp(-β₀ / β₁)

# Arguments
- `data::DataFrame`: Must have columns `:log_dose` (log of dose), `:n_dead`
  (number of deaths), and `:n_total` (total individuals tested).

# Returns
A `NamedTuple` with fields:
- `lc50::Float64`: Estimated LC50 on the original dose scale
- `slope::Float64`: Estimated slope coefficient (β₁)
"""
function fit_2pl(data::DataFrame)
    df = DataFrame(; log_dose = data.log_dose, proportion = data.n_dead ./ data.n_total)
    model = glm(
        @formula(proportion ~ log_dose),
        df,
        Binomial(),
        LogitLink();
        weights = FrequencyWeights(Float64.(data.n_total)),
    )
    β = coef(model)
    lc50 = exp(-β[1] / β[2])
    return (lc50 = lc50, slope = β[2])
end
