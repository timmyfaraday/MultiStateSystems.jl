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
        std = include(joinpath(_MSS.BASE_DIR,"test/elements/CAP_ac.jl"))
        
        solve!(std, MarkovProcess(), tsim = 1.0u"yr")
        prb_m = [np[end] for np in _MSS.get_sprop(std, :prob)]
        
        solve!(std, SemiMarkovProcess(), tsim = 1.0u"yr")  
        prb_s = [np[end] for np in _MSS.get_sprop(std, :prob)]
        
        @test isapprox(sum(prb_s), 1.0, rtol=1e-5)
        @test all(isapprox.(prb_m, prb_s, rtol=1e-5))
    end

    @testset "Optical Monitoring System" begin
        # Example considering a multi-state Markov model for an optical monitoring system
        # taken from:
        # > Mathematical formulation and numerical treatment based on transition frequency
        # densities and quadrature methods for non-homogeneous semi-Markov processes by 
        # M. das Chagas Moura, and E. L. Droguett (2009)

        stdᵒᵐˢ = include(joinpath(_MSS.BASE_DIR,"test/elements/optical_monitoring_system.jl"))
        solve!(stdᵒᵐˢ, SemiMarkovProcess(), tsim = 8760.0u"hr")
        prb_oms = [np[end] for np in _MSS.get_sprop(stdᵒᵐˢ, :prob)]
        @test isapprox(sum(prb_oms), 1.0, rtol=1e-4)
        @test isapprox(prb_oms[1], 0.8846075098143418, rtol = 1e-5)
        @test isapprox(prb_oms[2], 0.003709508099739578, rtol = 1e-5)
        @test isapprox(prb_oms[3], 0.001238170881240272, rtol = 1e-5)
        @test isapprox(prb_oms[4], 0.11045656347369365, rtol = 1e-5)
    end

    @testset "Semi-Analytical Example" begin
        # A semi-analytical example example of a multi-state semi-Markov model
        # taken from: 
        # > Mathematical formulation and numerical treatment based on transition frequency
        # densities and quadrature methods for non-homogeneous semi-Markov processes by 
        # M. das Chagas Moura, and E. L. Droguett (2009)
        stdˢᵃ = include(joinpath(_MSS.BASE_DIR,"test/elements/semi_analytical_example.jl"))
        solve!(stdˢᵃ, SemiMarkovProcess(), tsim = 4500.0u"hr", dt = 1.0u"hr")
        prb_sa = [np[end] for np in _MSS.get_sprop(stdˢᵃ, :prob)]

        # analytical solution
        dt      = 3.0u"hr"
        tsim    = 4500.0u"hr"
        t       = 0u"hr":dt:tsim
        dst1 = _MSS.Exponential(1000.0u"hr")
        dst2 = _MSS.Weibull(250.0u"hr", 1.5)
        PA1 = 1.0 .- _MSS.cdf.(dst1,t,t)
        PA2 = [quadgk(x -> _MSS.pdf.(dst1,x,x) * (1.0 - _MSS.cdf(dst2,nt.-x,nt)),zero(nt),nt,rtol=1e-8)[1] for nt in t]

        @test isapprox(sum(prb_sa), 1.0, rtol=1e-4)
        @test isapprox(prb_sa[1]+prb_sa[2], PA1[end]+PA2[end], rtol = 1e-5)
    end
end
