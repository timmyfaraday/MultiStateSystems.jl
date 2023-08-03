################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Semi Markov Process" begin

    @testset "AC/DC converter capacitance" begin
        # Example from: A semi-Markovian approach to evaluate the availability
        # of LVDC systems with integrated battery storage
        stdᶜᵃᵖ = include(joinpath(_MSS.BASE_DIR,"test/elements/CAP_ac.jl"))
        
        solve!(stdᶜᵃᵖ, MarkovProcess(), tsim = 1.0u"yr")
        prb_m = [np[end] for np in _MSS.get_sprop(stdᶜᵃᵖ, :prob)]
        
        solve!(stdᶜᵃᵖ, SemiMarkovProcess(), tsim = 1.0u"yr")  
        prb_s = [np[end] for np in _MSS.get_sprop(stdᶜᵃᵖ, :prob)]
        
        @test isapprox(sum(prb_m), 1.0, rtol=1e-6)
        @test all(isapprox.(prb_m, prb_s, rtol=1e-6))
    end

    @testset "Optical Monitoring System" begin
        # Example considering a multi-state Markov model for an optical monitoring system
        # taken from:
        # > A Multi-State Markov Model for a Short-Term Reliability Analysis of a Power
        # Generating Unit by A. Lisnianski, D. Elmakias, and H. Ben Haim (2012)

        stdᵒᵐˢ = include(joinpath(_MSS.BASE_DIR,"test/elements/optical_monitoring_system.jl"))
        solve!(stdᵒᵐˢ, SemiMarkovProcess(), tsim = 1.0u"yr")
        prb_oms = [np[end] for np in _MSS.get_sprop(stdᵒᵐˢ, :prob)]
        @test isapprox(sum(prb_oms), 1.0, rtol=1e-6)
    end

end
