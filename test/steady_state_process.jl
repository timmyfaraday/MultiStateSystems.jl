################################################################################
#  Copyright 2021, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Steady State Process" begin
    
    @testset "2003 System" begin
        # Example from: Reliability of Safety-Critical System - Theory and 
        # Applications (ch. 5)
        global λ = 5e-4u"1/hr"
        global μ = 1e-1u"1/hr"

        std²³ = include(joinpath(_MSS.BASE_DIR,"test/elements/2003_system.jl"))
        
        solve!(std²³, SteadyStateProcess())
        prb = [np[end] for np in _MSS.get_sprop(std²³, :prob)]
        sol = [ μ*(2*λ*λ + 2*λ*μ + μ*μ + λ*μ),
                3*λ*μ*(λ + μ),
                6*μ*λ*λ,
                6*λ*λ*λ] ./ (6*λ*λ*λ + 11*λ*λ*μ + 5*λ*μ*μ + μ*μ*μ + λ*μ*μ)

        @test isapprox(sum(prb), 1.0, rtol=1e-6)
        @test all(isapprox.(prb, sol, rtol=1e-6))
    end

    @testset "2003 System with Measurements Uncertainty" begin
        # Example from: Reliability of Safety-Critical System - Theory and 
        # Applications (ch. 5)
        global λ = (5e-4 ± 1e-5)u"1/hr"
        global μ = (1e-1 ± 1e-2)u"1/hr"

        std²³ = include(joinpath(_MSS.BASE_DIR,"test/elements/2003_system.jl"))

        solve!(std²³, SteadyStateProcess())
        prb = [np[end] for np in _MSS.get_sprop(std²³, :prob)]
        sol = [ μ*(2*λ*λ + 2*λ*μ + μ*μ + λ*μ),
                3*λ*μ*(λ + μ),
                6*μ*λ*λ,
                6*λ*λ*λ] ./ (6*λ*λ*λ + 11*λ*λ*μ + 5*λ*μ*μ + μ*μ*μ + λ*μ*μ)

        @test isapprox(sum(prb), 1.0, rtol=1e-6)
        @test all(isapprox.(prb, sol, rtol=1e-6))

        @test isapprox(_MSM.uncertainty(sum(prb)), 0.0, atol=1e-6)
    end
end