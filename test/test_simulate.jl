using StableRNGs: StableRNG
using DataFrames: DataFrame, nrow, names
using Statistics: mean

@testset "simulate_bioassay output structure" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)
    doses = 1:5
    df = simulate_bioassay(pop, doses, 3, 20; rng = StableRNG(1))

    @test df isa DataFrame
    @test Set(names(df)) == Set(["log_dose", "replicate", "n_dead", "n_total"])
    @test nrow(df) == length(doses) * 3
    @test all(==(20), df.n_total)
    @test all(x -> 0 <= x <= 20, df.n_dead)
    @test all(x -> x ∈ 1:3, df.replicate)
    @test all(>(0), df.n_total)
end

@testset "simulate_bioassay reproducibility" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)
    doses = 1:5

    df1 = simulate_bioassay(pop, doses, 2, 50; rng = StableRNG(42))
    df2 = simulate_bioassay(pop, doses, 2, 50; rng = StableRNG(42))
    @test df1.n_dead == df2.n_dead
    @test df1.log_dose == df2.log_dose

    df3 = simulate_bioassay(pop, doses, 2, 50; rng = StableRNG(99))
    @test df1.n_dead != df3.n_dead
end

@testset "simulate_bioassay sd_log_lc50 = 0" begin
    pop = UnimodalPopulation(log(5.0), 0.0, 10.0)
    df = simulate_bioassay(pop, [5.0], 1, 100; rng = StableRNG(7))
    @test df.n_dead[1] ≈ 50 atol = 15
end

@testset "simulate_bioassay extreme doses" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)

    # Near-zero mortality at very low dose
    df_low = simulate_bioassay(pop, [0.001], 5, 200; rng = StableRNG(10))
    @test mean(df_low.n_dead ./ df_low.n_total) < 0.05

    # Near-total mortality at very high dose
    df_high = simulate_bioassay(pop, [1000.0], 5, 200; rng = StableRNG(11))
    @test mean(df_high.n_dead ./ df_high.n_total) > 0.95
end

@testset "simulate_bioassay input validation" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)
    @test_throws ArgumentError simulate_bioassay(pop, [-1.0, 2.0], 2, 10)
    @test_throws ArgumentError simulate_bioassay(pop, [1.0, 2.0], 0, 10)
    @test_throws ArgumentError simulate_bioassay(pop, [1.0, 2.0], 2, 0)
end
