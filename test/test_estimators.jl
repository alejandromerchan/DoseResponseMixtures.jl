using StableRNGs: StableRNG
using DataFrames: names
using Statistics: mean

@testset "fit_2pl return type" begin
    pop = UnimodalPopulation(log(5.0), 0.0, 10.0)
    df = simulate_bioassay(pop, 1:10, 5, 100; rng = StableRNG(1))
    result = fit_2pl(df)
    @test result isa NamedTuple
    @test haskey(result, :lc50)
    @test haskey(result, :slope)
    @test result.lc50 isa Float64
    @test result.slope isa Float64
end

@testset "fit_2pl parameter recovery" begin
    true_lc50 = 5.0
    true_slope = 10.0
    pop = UnimodalPopulation(log(true_lc50), 0.0, true_slope)

    # Large dataset for tight estimates
    df = simulate_bioassay(pop, 1:10, 10, 1000; rng = StableRNG(123))
    result = fit_2pl(df)

    @test result.lc50 ≈ true_lc50 rtol = 0.05
    @test result.slope ≈ true_slope rtol = 0.10
end

@testset "fit_2pl does not mutate input" begin
    pop = UnimodalPopulation(log(5.0), 0.0, 10.0)
    df = simulate_bioassay(pop, 1:10, 3, 100; rng = StableRNG(7))
    cols_before = Set(names(df))
    fit_2pl(df)
    @test Set(names(df)) == cols_before
end
