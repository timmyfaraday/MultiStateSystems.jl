################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Universal Generating Operator" begin 

    @testset "Series-Parallel Network" begin
        # Example from: Multi-State System Reliability Analysis and Optimization
        # for Engineers and Industrial Managers (p. 11/209)
        val_sol = [0u"m^3/hr", 1500u"m^3/hr", 1800u"m^3/hr", 2000u"m^3/hr", 3500u"m^3/hr"]
        prb_sol = [0.01207138, 0.10333981359, 0.07371573369, 0.05304777100, 0.75782529972]

        global λ¹, μ¹ = 7.0u"1/yr", 100.0u"1/yr" 
        global λ², μ² = 10.0u"1/yr", 80.0u"1/yr"
        global λ³ᵃ, λ³ᵇ = 7.0u"1/yr", 10.0u"1/yr"
        global μ³ᵃ, μ³ᵇ = 120.0u"1/yr", 110.0u"1/yr"

        ntwᶠᵗˢ = include(joinpath(_MSS.BASE_DIR,"test/networks/flow_transmission_system.jl"))
        solve!(ntwᶠᵗˢ)
        val_ntw = ntwᶠᵗˢ.usr[1][:ugf].val
        prb_ntw = ntwᶠᵗˢ.usr[1][:ugf].prb

        @test isapprox.(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
    end

    @testset "Series-Parallel Network with Measurements Uncertainty" begin
        # Example from: Multi-State System Reliability Analysis and Optimization
        # for Engineers and Industrial Managers (p. 11/209)
        val_sol = [0u"m^3/hr", 1500u"m^3/hr", 1800u"m^3/hr", 2000u"m^3/hr", 3500u"m^3/hr"]
        prb_sol = [0.01207138, 0.10333981359, 0.07371573369, 0.05304777100, 0.75782529972]

        global λ¹, μ¹ = (7.0±1.0)u"1/yr", 100.0u"1/yr" 
        global λ², μ² = 10.0u"1/yr", (80.0±10.0)u"1/yr"
        global λ³ᵃ, λ³ᵇ = 7.0u"1/yr", (10.0±2.0)u"1/yr"
        global μ³ᵃ, μ³ᵇ = (120.0±5.0)u"1/yr", 110.0u"1/yr"

        ntwᶠᵗˢ = include(joinpath(_MSS.BASE_DIR,"test/networks/flow_transmission_system.jl"))
        solve!(ntwᶠᵗˢ)
        val_ntw = ntwᶠᵗˢ.usr[1][:ugf].val
        prb_ntw = ntwᶠᵗˢ.usr[1][:ugf].prb

        @test isapprox.(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))

        @test isapprox(_MSM.uncertainty(sum(prb_ntw)), 0.0, atol=1e-6)
    end

    @testset "Bridge Network" begin
        val_sol = [0.0u"m^3/hr", 1.0u"m^3/hr", 2.0u"m^3/hr"]
        prb_sol = [0.02152, 0.32238, 0.6561]

        ntwᵇⁿ = include(joinpath(_MSS.BASE_DIR,"test/networks/bridge_network.jl"))
        solve!(ntwᵇⁿ)
        val_ntw = ntwᵇⁿ.usr[1][:ugf].val
        prb_ntw = ntwᵇⁿ.usr[1][:ugf].prb

        @test isapprox.(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
    end

    @testset "Double Bridge Network" begin
        val_sol = [0.0u"m^3/hr", 1.0u"m^3/hr", 2.0u"m^3/hr"]
        prb_sol = [0.0102772, 0.1827198, 0.807003]

        ntwᵈᵇⁿ = include(joinpath(_MSS.BASE_DIR,"test/networks/double_bridge_network.jl"))
        solve!(ntwᵈᵇⁿ)
        val_ntw = ntwᵈᵇⁿ.usr[1][:ugf].val 
        prb_ntw = ntwᵈᵇⁿ.usr[1][:ugf].prb

        @test isapprox.(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
    end

    @testset "Disconnected Sources and Users" begin
        # Example from: Industrial Energy System Availability Management (p. 92)
        val_sol = [0.0u"MW", 1.0u"MW"]
        prb_sol = [0.1, 0.9]

        ntwᵈˢᶠ = include(joinpath(_MSS.BASE_DIR,"test/networks/double_single_feeder.jl"))
        solve!(ntwᵈˢᶠ)
        val_usr_1, prb_usr_1 = ntwᵈˢᶠ.usr[1][:ugf].val, ntwᵈˢᶠ.usr[1][:ugf].prb
        val_usr_2, prb_usr_2 = ntwᵈˢᶠ.usr[2][:ugf].val, ntwᵈˢᶠ.usr[2][:ugf].prb

        @test isapprox.(sum(prb_usr_1), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_usr_1, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_usr_1, rtol=1e-6))

        @test isapprox.(sum(prb_usr_2), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_usr_2, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_usr_2, rtol=1e-6))
    end

    @testset "Source Dependence" begin
        # Example from: Analytical Model for Availability Assessment of Large-
        # Scale Offshore Wind Farms including their Collector System (p. 5)
        val_sol = [0.0u"MW", 2.0u"MW", 4.0u"MW", 6.0u"MW", 8.0u"MW"]
        prb_sol = [0.307,0.0126,0.11907,0.10206,0.45927]

        ntwʷᶠ = include(joinpath(_MSS.BASE_DIR,"test/networks/wind_farm.jl"))
        solve!(ntwʷᶠ)
        val_ntw = ntwʷᶠ.usr[1][:ugf].val
        prb_ntw = ntwʷᶠ.usr[1][:ugf].prb

        @test isapprox.(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
    end

    @testset "Evaluation Dependence" begin
        ugfᵉ¹ = UGF(:power, [0.0u"MW",1.0u"MW",0.0u"MW",1.0u"MW"],
                            [0.1,0.2,0.3,0.4], rdc=false)
        ugfᵉ² = UGF(:power, [0.0u"MW",0.0u"MW",1.0u"MW",1.0u"MW"],
                            [0.1,0.2,0.3,0.4], rdc=false)

        # parallel src system: str = ugfᵉ¹.val + ugfᵉ².val 
        # 1) val = 0.0u"MW" → prb = 0.1 
        # 2) val = 1.0u"MW" → prb = 0.2 + 0.3 = 0.5
        # 2) val = 2.0u"MW" → prb = 0.4 
        val_sol = [0.0u"MW", 1.0u"MW", 2.0u"MW"] 
        prb_sol = [0.1, 0.5, 0.4]

        ntwˣ = Network()
        add_sources!(ntwˣ, node = [1,2],
                           ugf = [ugfᵉ¹,ugfᵉ²],
                           eval_dep = true)
        add_user!(ntwˣ, node = 3)
        add_components!(ntwˣ, edge = [(1,3),(2,3)])
        solve!(ntwˣ)
        val_ntw = ntwˣ.usr[1][:ugf].val
        prb_ntw = ntwˣ.usr[1][:ugf].prb

        @test isapprox(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))

        # series cmp system: str = min(ugfᵉ¹.val,ugfᵉ².val)
        # 1) val = 0.0u"MW" → prb = 0.1 + 0.2 + 0.3 = 0.6
        # 2) val = 1.0u"MW" → prb = 0.4 
        val_sol = [0.0u"MW", 1.0u"MW"]
        prb_sol = [0.6, 0.4]

        ntwˢ = Network()
        add_source!(ntwˢ, node = 1)
        add_user!(ntwˢ, node = 3)
        add_components!(ntwˢ, edge = [(1,2),(2,3)], 
                              ugf = [ugfᵉ¹,ugfᵉ²],
                              eval_dep = true)
        solve!(ntwˢ)
        val_ntw = ntwˢ.usr[1][:ugf].val
        prb_ntw = ntwˢ.usr[1][:ugf].prb

        @test isapprox(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
        
        # parallel cmp system: str = ugfᵉ¹.val + ugfᵉ².val
        # 1) val = 0.0u"MW" → prb = 0.1 
        # 2) val = 1.0u"MW" → prb = 0.2 + 0.3 = 0.5
        # 2) val = 2.0u"MW" → prb = 0.4 
        val_sol = [0.0u"MW", 1.0u"MW", 2.0u"MW"] 
        prb_sol = [0.1, 0.5, 0.4]

        ntwᵖ = Network()
        add_source!(ntwᵖ, node = 1)
        add_user!(ntwᵖ, node = 2)
        add_components!(ntwᵖ, edge = [(1,2),(1,2)], 
                              ugf = [ugfᵉ¹,ugfᵉ²],
                              eval_dep = true)
        solve!(ntwᵖ)
        val_ntw = ntwᵖ.usr[1][:ugf].val
        prb_ntw = ntwᵖ.usr[1][:ugf].prb

        @test isapprox(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
    end

end