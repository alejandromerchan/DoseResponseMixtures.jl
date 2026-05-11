using StableRNGs: StableRNG
using DataFrames: DataFrame, nrow, names
using Statistics: mean, std

@testset "HWEMixturePopulation constructor and validation" begin
    pop = HWEMixturePopulation(0.3, log(5.0), 2.0, 0.5, 0.2, 10.0)
    @test pop isa AbstractPopulation
    @test pop isa HWEMixturePopulation
    @test pop.allele_frequency === 0.3
    @test pop.mean_log_lc50_susceptible === log(5.0)
    @test pop.resistance_magnitude === 2.0
    @test pop.dominance === 0.5
    @test pop.sd_log_lc50 === 0.2
    @test pop.slope === 10.0

    # Integer args coerce to Float64
    pop2 = HWEMixturePopulation(0, 1, 2, 0, 0, 1)
    @test pop2.allele_frequency isa Float64

    # Boundary values are valid
    @test HWEMixturePopulation(0.0, 0.0, 0.0, 0.0, 0.0, 1.0) isa HWEMixturePopulation
    @test HWEMixturePopulation(1.0, 0.0, 0.0, 1.0, 0.0, 1.0) isa HWEMixturePopulation

    # allele_frequency out of range
    @test_throws ArgumentError HWEMixturePopulation(-0.01, 0.0, 1.0, 0.5, 0.2, 10.0)
    @test_throws ArgumentError HWEMixturePopulation(1.01, 0.0, 1.0, 0.5, 0.2, 10.0)

    # dominance out of range
    @test_throws ArgumentError HWEMixturePopulation(0.3, 0.0, 1.0, -0.01, 0.2, 10.0)
    @test_throws ArgumentError HWEMixturePopulation(0.3, 0.0, 1.0, 1.01, 0.2, 10.0)

    # sd_log_lc50 negative
    @test_throws ArgumentError HWEMixturePopulation(0.3, 0.0, 1.0, 0.5, -0.1, 10.0)

    # slope non-positive
    @test_throws ArgumentError HWEMixturePopulation(0.3, 0.0, 1.0, 0.5, 0.2, 0.0)
    @test_throws ArgumentError HWEMixturePopulation(0.3, 0.0, 1.0, 0.5, 0.2, -1.0)
end

@testset "genotype_means" begin
    pop = HWEMixturePopulation(0.3, 1.0, 2.0, 0.5, 0.2, 10.0)
    means = genotype_means(pop)
    @test means.SS === 1.0
    @test means.RR === 3.0
    @test means.SR === 2.0

    # h = 0: SR == SS
    pop0 = HWEMixturePopulation(0.5, 1.0, 2.0, 0.0, 0.2, 10.0)
    m0 = genotype_means(pop0)
    @test m0.SS == m0.SR
    @test m0.RR == 3.0

    # h = 1: SR == RR
    pop1 = HWEMixturePopulation(0.5, 1.0, 2.0, 1.0, 0.2, 10.0)
    m1 = genotype_means(pop1)
    @test m1.SR == m1.RR
    @test m1.SS == 1.0

    # h = 0.5: SR exactly halfway
    pop_add = HWEMixturePopulation(0.5, 1.0, 2.0, 0.5, 0.2, 10.0)
    m_add = genotype_means(pop_add)
    @test m_add.SR == (m_add.SS + m_add.RR) / 2

    # resistance_magnitude = 0: all means equal μ_SS
    pop_r0 = HWEMixturePopulation(0.5, 1.0, 0.0, 0.5, 0.2, 10.0)
    mr0 = genotype_means(pop_r0)
    @test mr0.SS == mr0.SR == mr0.RR
end

@testset "genotype_frequencies" begin
    pop = HWEMixturePopulation(0.3, 1.0, 2.0, 0.5, 0.2, 10.0)
    freqs = genotype_frequencies(pop)
    q = 0.3
    @test freqs.SS ≈ (1 - q)^2
    @test freqs.SR ≈ 2 * q * (1 - q)
    @test freqs.RR ≈ q^2
    @test freqs.SS + freqs.SR + freqs.RR ≈ 1.0

    # q = 0: only SS
    pop0 = HWEMixturePopulation(0.0, 1.0, 2.0, 0.5, 0.2, 10.0)
    f0 = genotype_frequencies(pop0)
    @test f0.SS == 1.0
    @test f0.SR == 0.0
    @test f0.RR == 0.0

    # q = 1: only RR
    pop1 = HWEMixturePopulation(1.0, 1.0, 2.0, 0.5, 0.2, 10.0)
    f1 = genotype_frequencies(pop1)
    @test f1.SS == 0.0
    @test f1.SR == 0.0
    @test f1.RR == 1.0
end

@testset "limit cases" begin
    n = 10_000
    μ_SS = log(5.0)
    Δ = 2.0
    σ = 0.3

    # q = 0: all individuals should be :SS
    pop0 = HWEMixturePopulation(0.0, μ_SS, Δ, 0.5, σ, 10.0)
    genotypes0, lts0 =
        DoseResponseMixtures._sample_genotypes_and_log_tolerances(pop0, n, StableRNG(1))
    @test all(==(:SS), genotypes0)
    @test mean(lts0) ≈ μ_SS atol = 0.05

    # q = 1: all individuals should be :RR
    pop1 = HWEMixturePopulation(1.0, μ_SS, Δ, 0.5, σ, 10.0)
    genotypes1, lts1 =
        DoseResponseMixtures._sample_genotypes_and_log_tolerances(pop1, n, StableRNG(2))
    @test all(==(:RR), genotypes1)
    @test mean(lts1) ≈ μ_SS + Δ atol = 0.05

    # resistance_magnitude = 0: all means equal μ_SS, samples ~Normal(μ_SS, σ)
    pop_r0 = HWEMixturePopulation(0.5, μ_SS, 0.0, 0.5, σ, 10.0)
    lts_r0 = sample_log_tolerances(pop_r0, n; rng = StableRNG(3))
    @test mean(lts_r0) ≈ μ_SS atol = 0.05
    @test std(lts_r0) ≈ σ atol = 0.05
end

@testset "HWE verification" begin
    q = 0.3
    n = 100_000
    pop = HWEMixturePopulation(q, log(5.0), 2.0, 0.5, 0.2, 10.0)
    genotypes, _ =
        DoseResponseMixtures._sample_genotypes_and_log_tolerances(pop, n, StableRNG(42))
    obs_SS = count(==(:SS), genotypes) / n
    obs_SR = count(==(:SR), genotypes) / n
    obs_RR = count(==(:RR), genotypes) / n

    @test obs_SS ≈ (1 - q)^2 atol = 0.01
    @test obs_SR ≈ 2 * q * (1 - q) atol = 0.01
    @test obs_RR ≈ q^2 atol = 0.01
end

@testset "dominance behavior" begin
    # h = 0: SS and SR have same mean
    pop_h0 = HWEMixturePopulation(0.5, 1.0, 2.0, 0.0, 0.2, 10.0)
    m_h0 = genotype_means(pop_h0)
    @test m_h0.SS == m_h0.SR

    # h = 1: SR and RR have same mean
    pop_h1 = HWEMixturePopulation(0.5, 1.0, 2.0, 1.0, 0.2, 10.0)
    m_h1 = genotype_means(pop_h1)
    @test m_h1.SR == m_h1.RR

    # h = 0.5: SR exactly halfway between SS and RR
    pop_h5 = HWEMixturePopulation(0.5, 1.0, 2.0, 0.5, 0.2, 10.0)
    m_h5 = genotype_means(pop_h5)
    @test m_h5.SR == (m_h5.SS + m_h5.RR) / 2
end

@testset "HWEMixturePopulation sample_log_tolerances" begin
    pop = HWEMixturePopulation(0.3, log(5.0), 2.0, 0.5, 0.2, 10.0)

    s = sample_log_tolerances(pop, 100; rng = StableRNG(42))
    @test s isa Vector{Float64}
    @test length(s) == 100

    # Reproducibility
    s1 = sample_log_tolerances(pop, 200; rng = StableRNG(7))
    s2 = sample_log_tolerances(pop, 200; rng = StableRNG(7))
    @test s1 == s2

    # Different seeds produce different results
    s3 = sample_log_tolerances(pop, 200; rng = StableRNG(99))
    @test s1 != s3

    # Input validation
    @test_throws ArgumentError sample_log_tolerances(pop, 0)
    @test_throws ArgumentError sample_log_tolerances(pop, -1)
end

@testset "HWEMixturePopulation simulate_bioassay" begin
    pop = HWEMixturePopulation(0.3, log(5.0), 2.0, 0.5, 0.2, 10.0)
    doses = 1:5

    df = simulate_bioassay(pop, doses, 3, 20; rng = StableRNG(1))

    @test df isa DataFrame
    @test Set(names(df)) == Set(["log_dose", "replicate", "n_dead", "n_total"])
    @test nrow(df) == length(doses) * 3
    @test all(==(20), df.n_total)
    @test all(x -> 0 <= x <= 20, df.n_dead)

    # No genotype column
    @test !("genotype" in names(df))

    # Reproducibility
    df1 = simulate_bioassay(pop, doses, 2, 50; rng = StableRNG(42))
    df2 = simulate_bioassay(pop, doses, 2, 50; rng = StableRNG(42))
    @test df1.n_dead == df2.n_dead
    @test df1.log_dose == df2.log_dose

    df3 = simulate_bioassay(pop, doses, 2, 50; rng = StableRNG(99))
    @test df1.n_dead != df3.n_dead

    # Input validation
    @test_throws ArgumentError simulate_bioassay(pop, [-1.0, 2.0], 2, 10)
    @test_throws ArgumentError simulate_bioassay(pop, [1.0, 2.0], 0, 10)
    @test_throws ArgumentError simulate_bioassay(pop, [1.0, 2.0], 2, 0)
end
