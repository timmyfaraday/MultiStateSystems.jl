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
        # for Engineers and Industrial Managers (p. 11)
        val_sol = Any[0u"m^3/hr", 1500u"m^3/hr", 1800u"m^3/hr", 2000u"m^3/hr", 3500u"m^3/hr"]
        prb_sol = Any[0.172, 0.288, 0.12, 0.084, 0.336]

        ntwᶠᵗˢ = include(joinpath(_MSS.BASE_DIR,"test/networks/flow_transmission_system.jl"))
        solve!(ntwᶠᵗˢ)
        val_ntw = ntwᶠᵗˢ.usr[1][:ugf].val
        prb_ntw = ntwᶠᵗˢ.usr[1][:ugf].prb

        @test isapprox.(sum(prb_ntw), 1.0, rtol=1e-6)
        @test all(isapprox.(val_sol, val_ntw, rtol=1e-6))
        @test all(isapprox.(prb_sol, prb_ntw, rtol=1e-6))
    end

    @testset "Disconnected Sources and Users" begin
        # Example from: Industrial Energy System Availability Management (p. 92)
        val_sol = Any[0.0u"MW", 1.0u"MW"]
        prb_sol = Any[0.1, 0.9]

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
        val_sol = Any[0.0u"MW", 2.0u"MW", 4.0u"MW", 6.0u"MW", 8.0u"MW"]
        prb_sol = Any[0.307,0.0126,0.11907,0.10206,0.45927]

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