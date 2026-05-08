using Random: AbstractRNG, default_rng
using Distributions: Normal
using DataFrames: DataFrame

"""
    simulate_bioassay(pop, doses, n_replicates, n_per_replicate; rng)

Simulate a dose-response bioassay for a `UnimodalPopulation`.

For each (dose, replicate) combination, draws `n_per_replicate` individual
log-LC50 values from the population tolerance distribution, then samples
Bernoulli mortality outcomes using the individual-level log-logistic model:

    p = 1 / (1 + exp(-slope * (log(dose) - log_lc50_individual)))

# Arguments
- `pop::UnimodalPopulation`: Population tolerance distribution
- `doses`: Iterable of positive dose values
- `n_replicates::Integer`: Number of replicates per dose level
- `n_per_replicate::Integer`: Number of individuals per replicate
- `rng::AbstractRNG`: Random number generator (default: `Random.default_rng()`)

# Returns
A `DataFrame` with columns `:log_dose`, `:replicate`, `:n_dead`, `:n_total`.
"""
function simulate_bioassay(
    pop::UnimodalPopulation,
    doses,
    n_replicates::Integer,
    n_per_replicate::Integer;
    rng::AbstractRNG = default_rng(),
)
    all(>(0), doses) || throw(ArgumentError("all doses must be positive"))
    n_replicates > 0 || throw(ArgumentError("n_replicates must be positive"))
    n_per_replicate > 0 || throw(ArgumentError("n_per_replicate must be positive"))

    dose_vec = collect(doses)
    n_rows = length(dose_vec) * n_replicates

    log_dose_col = Vector{Float64}(undef, n_rows)
    replicate_col = Vector{Int}(undef, n_rows)
    n_dead_col = Vector{Int}(undef, n_rows)
    n_total_col = fill(Int(n_per_replicate), n_rows)

    row = 0
    for dose in dose_vec
        log_dose = log(Float64(dose))
        for rep = 1:n_replicates
            row += 1
            # Draw individual log-LC50 values (degenerate case handled explicitly
            # because Normal(μ, 0) is not a valid distribution)
            log_lc50s = if pop.sd_log_lc50 > 0
                rand(rng, Normal(pop.mean_log_lc50, pop.sd_log_lc50), n_per_replicate)
            else
                fill(pop.mean_log_lc50, n_per_replicate)
            end
            n_dead = 0
            for log_lc50 in log_lc50s
                p = 1 / (1 + exp(-pop.slope * (log_dose - log_lc50)))
                n_dead += rand(rng) < p
            end
            log_dose_col[row] = log_dose
            replicate_col[row] = rep
            n_dead_col[row] = n_dead
        end
    end

    return DataFrame(;
        log_dose = log_dose_col,
        replicate = replicate_col,
        n_dead = n_dead_col,
        n_total = n_total_col,
    )
end
