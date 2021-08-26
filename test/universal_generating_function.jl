################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Universal Generating Function" begin 

    @testset "Optional Reduction" begin
        ugfᵍᵉⁿ = UGF(:power, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7])
        @test isapprox(ugfᵍᵉⁿ.val, [0.0u"MW",2.0u"MW"], rtol=1e-6)
        @test isapprox(ugfᵍᵉⁿ.prb, [0.3,0.7], rtol=1e-6)

        ugfᵍᵉⁿ = UGF(:power, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7], rdc=false)
        @test isapprox(ugfᵍᵉⁿ.val, [0.0u"MW",0.0u"MW",2.0u"MW"], rtol=1e-6)
        @test isapprox(ugfᵍᵉⁿ.prb, [0.1, 0.2, 0.7], rtol=1e-6)
    end

    @testset "STD → UGF" begin
        stdᵍᵉⁿ = solvedSTD(prob  = [0.1,0.2,0.7],
                           power = [0.0u"MW",0.0u"MW",2.0u"MW"])
        ugfᵍᵉⁿ = UGF(:power, stdᵍᵉⁿ)
        @test isapprox(ugfᵍᵉⁿ.val, [0.0u"MW",2.0u"MW"], rtol=1e-6)
        @test isapprox(ugfᵍᵉⁿ.prb, [0.3,0.7], rtol=1e-6)
    end
end