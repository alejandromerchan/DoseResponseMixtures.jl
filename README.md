# DoseResponseMixtures

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://alejandromerchan.github.io/DoseResponseMixtures.jl/stable)
[![Development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://alejandromerchan.github.io/DoseResponseMixtures.jl/dev)
[![Test workflow status](https://github.com/alejandromerchan/DoseResponseMixtures.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/alejandromerchan/DoseResponseMixtures.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/alejandromerchan/DoseResponseMixtures.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/alejandromerchan/DoseResponseMixtures.jl)
[![Lint workflow Status](https://github.com/alejandromerchan/DoseResponseMixtures.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/alejandromerchan/DoseResponseMixtures.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/alejandromerchan/DoseResponseMixtures.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/alejandromerchan/DoseResponseMixtures.jl/actions/workflows/Docs.yml?query=branch%3Amain)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

> **Status:** Early development. The package is being built incrementally as a side project; APIs will change without warning until a 1.0 release.

## Overview

`DoseResponseMixtures.jl` is a Julia package for analyzing dose-response bioassays under the assumption that tolerance to a stressor (typically a pesticide) is a quantitative trait whose distribution in a population may be a *mixture* of components rather than a single unimodal distribution.

The classical bioassay infrastructure — probit analysis, the 2-parameter log-logistic model, and similar tools — assumes a single tolerance distribution. This assumption fails precisely when it matters most: during the spread of a resistance allele through a population, when the population is structurally a mixture of susceptible, heterozygous, and resistant individuals. The same failure occurs by design in refuge-based resistance management strategies (e.g., for Bt-resistant *Helicoverpa zea*), where the population is *intended* to be a mixture of genotypes.

This package provides:

- A **mechanistic generative model** of bioassay data, in which a population is parameterized by allele frequency, genotype-specific mean tolerances, dominance, and within-genotype variance, and bioassay outcomes follow from individual-level Bernoulli mortality.
- A **suite of estimators** ranging from the classical 2-parameter log-logistic to mixture models with explicit quantitative-genetic structure under Hardy-Weinberg equilibrium.
- **Evaluation tooling** for comparing estimators on bias, precision, and detection power across realistic experimental designs.

## Scientific motivation

A more detailed discussion of the modeling choices is forthcoming in the documentation. In brief:

1. Pesticide tolerance has a polygenic basis and is reasonably modeled, on a log scale, as approximately normally distributed within a genotype class.
2. During the spread of a major-effect resistance allele, or in any population maintained as a genotype mixture by management practice, the phenotypic tolerance distribution is *multimodal*.
3. Estimators built on a unimodal assumption produce systematically biased point estimates and slope estimates under these conditions, and standard goodness-of-fit checks frequently fail to detect the problem at realistic sample sizes.
4. Phenomenological extensions of the classical model (e.g., the 5-parameter log-logistic) improve fit but do not recover the underlying genetic structure.
5. Explicit mixture models, particularly those constrained by Hardy-Weinberg expectations, recover biologically meaningful parameters (allele frequency, dominance, resistance magnitude) and provide a principled framework for resistance monitoring.

## Installation

The package is not yet registered. To install the development version:

```julia
using Pkg
Pkg.add(url = "https://github.com/alejandromerchan/DoseResponseMixtures.jl")
```

## Roadmap

- [ ] Generative model: single-population (unimodal) baseline matching the original R simulations
- [ ] Generative model: HWE-structured three-genotype population with dominance
- [ ] Estimator: classical 2-parameter log-logistic via GLM
- [ ] Estimator: 2-parameter log-logistic with background mortality
- [ ] Estimator: agnostic two-component mixture (Turing.jl)
- [ ] Estimator: HWE-constrained three-component mixture with dominance (Turing.jl)
- [ ] Grid-runner for factorial simulation experiments
- [ ] Evaluation utilities: bias, RMSE, coverage, detection power
- [ ] Plotting recipes for standard diagnostic figures
- [ ] Documentation and worked examples
- [ ] Manuscript reproduction script

## License

MIT. See [LICENSE](LICENSE).
