using StableRNGs: StableRNG
using Statistics: mean, std

@testset "sample_log_tolerances output structure" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)
    s = sample_log_tolerances(pop, 100; rng = StableRNG(42))

    @test s isa Vector{Float64}
    @test length(s) == 100
end

@testset "sample_log_tolerances reproducibility" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)

    s1 = sample_log_tolerances(pop, 50; rng = StableRNG(42))
    s2 = sample_log_tolerances(pop, 50; rng = StableRNG(42))
    @test s1 == s2

    s3 = sample_log_tolerances(pop, 50; rng = StableRNG(99))
    @test s1 != s3
end

@testset "sample_log_tolerances statistical properties" begin
    pop = UnimodalPopulation(log(5.0), 0.3, 10.0)
    s = sample_log_tolerances(pop, 10_000; rng = StableRNG(1))

    @test mean(s) ≈ pop.mean_log_lc50 atol = 0.05
    @test std(s) ≈ pop.sd_log_lc50 atol = 0.05
end

@testset "sample_log_tolerances sd_log_lc50 = 0" begin
    pop = UnimodalPopulation(log(5.0), 0.0, 10.0)
    s = sample_log_tolerances(pop, 50; rng = StableRNG(1))

    @test all(==(pop.mean_log_lc50), s)
end

@testset "sample_log_tolerances input validation" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)
    @test_throws ArgumentError sample_log_tolerances(pop, 0)
    @test_throws ArgumentError sample_log_tolerances(pop, -1)
end
