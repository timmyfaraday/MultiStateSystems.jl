using MultiStateSystems
using Test

@testset "State Transition Diagram" begin
    stdᵗᵉˢᵗ = STD()
    @test add_state!(stdᵗᵉˢᵗ) == true
end
