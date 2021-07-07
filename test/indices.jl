################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Indices" begin
    
    @testset "EENS" begin
        # Example from: Analytical Model for Availability Assessment of Large-
        # Scale Offshore Wind Farms including their Collector System (p. 6)
        ugf = UGF(:power,   [0.0u"MW", 2.0u"MW", 4.0u"MW", 6.0u"MW", 8.0u"MW"], 
                            [0.307000, 0.012600, 0.119070, 0.102060, 0.459270])
        pcc = _MSS.PropDict(:ugf => ugf)
        @test isapprox(_MSS.EENS(pcc), 28137.12u"MWh")
    end

    @testset "GRA" begin
        # Example from: Analytical Model for Availability Assessment of Large-
        # Scale Offshore Wind Farms including their Collector System (p. 6)
        ugf = UGF(:power,   [0.0u"MW", 1.0u"MW", 2.0u"MW", 3.0u"MW", 4.0u"MW"], 
                            [0.010000, 0.018000, 0.170100, 0.145800, 0.656100])
        ntw = _MSS.PropDict(:ugf => ugf)
        @test isapprox(_MSS.GRO(ntw, 0.7), 2.8u"MW")
        @test isapprox(_MSS.GRA(ntw, 0.7), 0.8019)
    end

end