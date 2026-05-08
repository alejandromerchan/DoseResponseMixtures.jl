@testset "UnimodalPopulation constructors" begin
    pop = UnimodalPopulation(log(5.0), 0.2, 10.0)
    @test pop.mean_log_lc50 ≈ log(5.0)
    @test pop.sd_log_lc50 ≈ 0.2
    @test pop.slope ≈ 10.0
    @test pop isa AbstractPopulation
    @test pop isa UnimodalPopulation
end

@testset "UnimodalPopulation field types" begin
    pop = UnimodalPopulation(1, 0, 2)
    @test pop.mean_log_lc50 isa Float64
    @test pop.sd_log_lc50 isa Float64
    @test pop.slope isa Float64
end

@testset "UnimodalPopulation validation" begin
    @test_throws ArgumentError UnimodalPopulation(0.0, -0.1, 10.0)
    @test_throws ArgumentError UnimodalPopulation(0.0, 0.2, 0.0)
    @test_throws ArgumentError UnimodalPopulation(0.0, 0.2, -1.0)

    # boundary: sd == 0 is valid (homogeneous population)
    @test UnimodalPopulation(1.0, 0.0, 5.0) isa UnimodalPopulation
end
