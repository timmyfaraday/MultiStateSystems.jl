################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Markov Process" begin

    @testset "Coal Fired Power Generating Unit" begin
        # Example from: A Multi-State Markov Model for a Short-Term Reliability
        # Analysis of a Power Generating Unit (p. 4)
        stdᶜᶠ = include(joinpath(_MSS.BASE_DIR,"test/elements/coal_fired_power_generating_unit.jl"))
        
        solve!(stdᶜᶠ, MarkovProcess(), tsim = 1.0u"yr")
        prb_m = [np[end] for np in _MSS.get_sprop(stdᶜᶠ, :prob)]
        
        solve!(stdᶜᶠ, SteadyStateProcess())  
        prb_s = [np[end] for np in _MSS.get_sprop(stdᶜᶠ, :prob)]
        
        @test isapprox(sum(prb_m), 1.0, rtol=1e-6)
        @test all(isapprox.(prb_m, prb_s, rtol=1e-6))
    end

end