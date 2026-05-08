"""
    DoseResponseMixtures

Dose-response bioassay analysis under mixture-distributed tolerance, with a
focus on populations whose phenotypic tolerance distribution arises from
quantitative-genetic mixtures of genotype classes (e.g., during the spread
of a resistance allele or under refuge-based resistance management).

This package is under active development; the public API is unstable.
"""
module DoseResponseMixtures

# Public API will be exported here as it is developed.
# Initial development order (see README roadmap):
#   1. Generative model (single-population baseline)
#   2. Generative model (HWE three-genotype mixture)
#   3. Classical 2PL estimator
#   4. 2PL with background mortality
#   5. Two-component mixture estimator (Turing.jl)
#   6. HWE three-component mixture estimator (Turing.jl)

"""
    hello_world()

Placeholder function from the package template. Will be removed once the
generative model is implemented.
"""
function hello_world()
    return "Hello, World!"
end

end
