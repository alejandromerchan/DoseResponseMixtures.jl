using Test
using DoseResponseMixtures

@testset "DoseResponseMixtures.jl" begin
    include("test_population.jl")
    include("test_simulate.jl")
    include("test_estimators.jl")
    include("test_sample.jl")
end
