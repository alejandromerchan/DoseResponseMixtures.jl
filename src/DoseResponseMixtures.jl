"""
    DoseResponseMixtures

Dose-response bioassay analysis under mixture-distributed tolerance, with a
focus on populations whose phenotypic tolerance distribution arises from
quantitative-genetic mixtures of genotype classes (e.g., during the spread
of a resistance allele or under refuge-based resistance management).

This package is under active development; the public API is unstable.
"""
module DoseResponseMixtures

include("population.jl")
include("simulate.jl")
include("estimators.jl")

export AbstractPopulation, UnimodalPopulation, HWEMixturePopulation
export genotype_means, genotype_frequencies
export simulate_bioassay, sample_log_tolerances
export fit_2pl

end
