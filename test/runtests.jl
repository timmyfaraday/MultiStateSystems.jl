using MultiStateSystems
using Test

@testset "State Transition Diagram" begin
    stdᵗᵉˢᵗ = STD()
    @test length(stdᵗᵉˢᵗ.sprops) == 0
end
