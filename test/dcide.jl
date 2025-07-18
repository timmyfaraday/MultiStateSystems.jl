################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using Test
using Unitful
using MultiStateSystems

@testset "DCIDE Module" begin
    dcide_bus_types = ["bus", "panel"]
    dcide_cmp_types = ["panel", "switch", "protection", "converter", "transformer", "bus"]
    dcide_src_types = ["utility", "battery", "PV"]
    dcide_usr_types = ["motor"]

    dcide_cmp_keys = [:id, :name, :std, :type, :node, :edge, :bidirectional]
    dcide_src_keys = [:id, :name, :std, :type, :node]
    dcide_usr_keys = [:id, :name, :type, :node]

    @testset "Filter Functions" begin
        # Create test component data
        

        test_cmp = (
            id = [1, 2, 3, 4],
            name = ["comp1", "comp2", "comp3", "comp4"],
            std = [nothing, nothing, nothing, nothing],
            type = ["bus", "panel", "switch", "cable"],
            node = [1, 2, nothing, nothing],
            edge = [nothing, nothing, (1, 2), (2, 3)],
            bidirectional = [false, false, true, false]
        )
        
        @testset "filter_nodal_cmp" begin
            nodal_cmp = _MSS.filter_nodal_cmp(test_cmp)
            @test length(nodal_cmp.id) == 2  # Only components with non-nothing nodes
            @test nodal_cmp.id == [1, 2]
            @test nodal_cmp.name == ["comp1", "comp2"]
            @test !haskey(nodal_cmp, :edge)  # edge key should be excluded
        end
        
        @testset "filter_bidirectional_edge_cmp" begin
            bidir_cmp = _MSS.filter_bidirectional_edge_cmp(test_cmp)
            @test length(bidir_cmp.id) == 1  # Only bidirectional edge components
            @test bidir_cmp.id == [3]
            @test bidir_cmp.edge == [(1, 2)]
            @test !haskey(bidir_cmp, :node)  # node key should be excluded
        end
        
        @testset "filter_unidirectional_edge_cmp" begin
            unidir_cmp = _MSS.filter_unidirectional_edge_cmp(test_cmp)
            @test length(unidir_cmp.id) == 1  # Only unidirectional edge components
            @test unidir_cmp.id == [4]
            @test unidir_cmp.edge == [(2, 3)]
            @test !haskey(unidir_cmp, :node)  # node key should be excluded
        end
        
        @testset "filter_panel_cmp" begin
            panel_cmp = _MSS.filter_panel_cmp(test_cmp)
            @test length(panel_cmp.id) == 2  # Only bus/panel types with non-nothing nodes
            @test panel_cmp.id == [1, 2]
            @test all(t -> t in dcide_bus_types, panel_cmp.type)
            @test !haskey(panel_cmp, :edge)  # edge key should be excluded
        end
    end
    
    @testset "get_std Function" begin
        # Test case 1: Component without specifications
        component_no_specs = Dict("id" => "test", "componentType" => "switch")
        std1 = _MSS.get_std(component_no_specs)
        
        @test std1 isa MultiStateSystems.AbstractSTD
        @test length(std1.sprops[1][:prob]) == 1
        @test std1.sprops[1][:prob] == 1
        @test length(std1.sprops[1][:power]) == 1
        @test std1.sprops[1][:power] == (Inf)u"MW"
        
        # Test case 2: Component with specifications but no availability
        component_no_av = Dict(
            "id" => "test",
            "componentType" => "switch",
            "specifications" => Dict("other" => "data")
        )
        std2 = _MSS.get_std(component_no_av)
        
        @test std2 isa MultiStateSystems.AbstractSTD
        @test std2.sprops[1][:prob] == 1
        @test std2.sprops[1][:power] == (Inf)u"MW"
        
        # Test case 3: Component with availability value A
        component_with_A = Dict(
            "id" => "test",
            "componentType" => "switch",
            "specifications" => Dict(
                "availability" => Dict(:A => 0.95),
                "electrical" => Dict(
                    "ports" => [Dict(
                        "AC" => Dict("power" => Dict("nom" => 1000))
                    )]
                )
            )
        )
        std3 = _MSS.get_std(component_with_A)
        
        @test std3 isa MultiStateSystems.AbstractSTD
        @test length(std3.sprops) == 2
        @test std3.sprops[1][:prob] ≈ 0.95
        @test std3.sprops[2][:prob] ≈ 0.05
        @test std3.sprops[1][:power] == 1000u"W"
        @test std3.sprops[2][:power] == 0.0u"W"
        
        # Test case 4: Component with DC power specification
        component_with_DC = Dict(
            "id" => "test",
            "componentType" => "converter",
            "specifications" => Dict(
                "availability" => Dict(:A => 0.99),
                "electrical" => Dict(
                    "ports" => [Dict(
                        "DC" => Dict("power" => Dict("nom" => 2000))
                    )]
                )
            )
        )
        std4 = _MSS.get_std(component_with_DC)
        
        @test std4 isa MultiStateSystems.AbstractSTD
        @test std4.sprops[1][:power] == 2000u"W"
        @test std4.sprops[2][:power] == 0.0u"W"
        
        # Test case 5: Component with failure/repair distributions
        component_with_distr = Dict(
            "id" => "test",
            "componentType" => "switch",
            "specifications" => Dict(
                "availability" => Dict(
                    :failure => Dict(:λ => 0.001, :k => nothing, :distr => _MSS.Exponential),
                    :repair => Dict(:μ => 0.1, :σ => nothing, :distr => _MSS.Exponential)
                )
            )
        )
        std5 = _MSS.get_std(component_with_distr)
        
        @test std5 isa MultiStateSystems.AbstractSTD
        @test _MSS.ns(std5) == 2  # Should have 2 states: Available and Unavailable
        @test _MSS.nt(std5) == 2  # Should have 2 transitions: (1,2) and (2,1)
    end
    
    @testset "init_elements Function" begin
        # Create test JSON data
        test_json = Dict(
            "componentInstances" => [
                Dict(
                    "id" => "PV1",
                    "componentType" => "PV",
                    "configuration" => Dict("ports" => ["port1"])
                ),
                Dict(
                    "id" => "Motor1",
                    "componentType" => "motor",
                    "configuration" => Dict("ports" => ["port1"])
                ),
                Dict(
                    "id" => "Switch1",
                    "componentType" => "switch",
                    "configuration" => Dict("ports" => ["port1", "port2"])
                ),
                Dict(
                    "id" => "Bus1",
                    "componentType" => "bus",
                    "configuration" => Dict("ports" => ["port1", "port2", "port3"])
                )
            ],
            "connections" => [
                Dict(
                    "id" => "Cable1",
                    "from" => Dict("componentInstanceId" => "PV1"),
                    "to" => Dict("componentInstanceId" => "Bus1")
                )
            ]
        )
        
        cmp, src, usr = _MSS.init_elements(test_json)
        
        @testset "Source Elements" begin
            @test length(src.id) == 1
            @test src.name[1] == "PV1"
            @test src.type[1] == "PV"
            @test src.node[1] == 0
        end
        
        @testset "User Elements" begin
            @test length(usr.id) == 1
            @test usr.name[1] == "Motor1"
            @test usr.type[1] == "motor"
            @test usr.node[1] == 0
        end
        
        @testset "Component Elements" begin
            # Should have Switch1, Bus1, and Cable1
            @test length(cmp.id) == 3
            @test "Switch1" in cmp.name
            @test "Bus1" in cmp.name
            @test "Cable1" in cmp.name
            
            # Check edge/node assignment logic
            switch_idx = findfirst(x -> x == "Switch1", cmp.name)
            bus_idx = findfirst(x -> x == "Bus1", cmp.name)
            cable_idx = findfirst(x -> x == "Cable1", cmp.name)
            
            @test cmp.edge[switch_idx] == (0, 0)  # 2 ports -> edge component
            @test cmp.node[switch_idx] === nothing
            
            @test cmp.node[bus_idx] == 0  # >2 ports -> nodal component
            @test cmp.edge[bus_idx] === nothing
            
            @test cmp.edge[cable_idx] == (0, 0)  # connection -> edge component
            @test cmp.node[cable_idx] === nothing
        end
    end
    
    @testset "connect_elements! Function" begin
        # Create test data
        cmp = (
            id = [1, 2],
            name = ["Cable1", "Bus1"],
            std = [nothing, nothing],
            type = ["cable", "bus"],
            node = [nothing, 0],
            edge = [(0, 0), nothing],
            bidirectional = [true, false]
        )
        
        src = (
            id = [1],
            name = ["PV1"],
            std = [nothing],
            type = ["PV"],
            node = [0]
        )
        
        usr = (
            id = [1],
            name = ["Motor1"],
            type = ["motor"],
            node = [0]
        )
        
        test_json = Dict(
            "connections" => [
                Dict(
                    "id" => "Cable1",
                    "from" => Dict("componentInstanceId" => "PV1"),
                    "to" => Dict("componentInstanceId" => "Bus1")
                )
            ]
        )
        
        # Convert to mutable arrays for the function to modify
        cmp_mut = (
            id = copy(cmp.id),
            name = copy(cmp.name),
            std = copy(cmp.std),
            type = copy(cmp.type),
            node = copy(cmp.node),
            edge = copy(cmp.edge),
            bidirectional = copy(cmp.bidirectional)
        )
        
        src_mut = (
            id = copy(src.id),
            name = copy(src.name),
            std = copy(src.std),
            type = copy(src.type),
            node = copy(src.node)
        )
        
        usr_mut = (
            id = copy(usr.id),
            name = copy(usr.name),
            type = copy(usr.type),
            node = copy(usr.node)
        )
        
        _MSS.connect_elements!(cmp_mut, src_mut, usr_mut, test_json)
        
        # After connection, nodes should be assigned non-zero values
        @test src_mut.node[1] != 0
        @test cmp_mut.node[2] != 0  # Bus1 should get a node assignment
        @test cmp_mut.edge[1] != (0, 0)  # Cable1 should get proper edge assignment
    end
    
    @testset "build_network Function" begin
        # Create simple test data
        cmp = (
            id = [1],
            name = ["Bus1"],
            std = [_MSS.get_std(Dict())],
            type = ["bus"],
            node = [1],
            edge = [nothing],
            bidirectional = [false]
        )
        
        src = (
            id = [1],
            name = ["PV1"],
            std = [_MSS.get_std(Dict())],
            type = ["PV"],
            node = [2]
        )
        
        usr = (
            id = [1],
            name = ["Motor1"],
            type = ["motor"],
            node = [3]
        )
        
        ntw = _MSS.build_network(cmp, src, usr)
        
        @test ntw isa MultiStateSystems.AbstractNetwork
        @test length(ntw.src) >= 1  # Should have at least the source
        @test length(ntw.usr) >= 1  # Should have at least the user
    end
    
    @testset "extend_json! Function" begin
        # Create test network and JSON
        test_json = Dict{String,Any}(
            "componentInstances" => [
            Dict{String,Any}("id" => 1, "name" => "PV1"),
            Dict{String,Any}("id" => 2, "name" => "Bus1", "componentType" => "bus")
            ],
            "connections" => [
            Dict{String,Any}("id" => 1, "name" => "Cable1")
            ]
        );
        
        # Create mock network with minimal structure
        ntw = Network()

        add_source!(ntw, node = 1, std = _MSS.get_std(test_json["componentInstances"][1]), id = 1)
        add_component!(ntw, edge = (1,2), std = _MSS.get_std(test_json["connections"][1]), type = "cable", id = 1)
        add_user!(ntw, node = 2, id = 1)

        solve!(ntw);

        # This should not throw an error
        @test_nowarn _MSS.extend_json!(test_json, ntw)
    end
    
    @testset "solve! Function Integration" begin
        # Create minimal test JSON with explicit Dict{String,Any} type
        test_json = Dict{String,Any}(
            "componentInstances" => [
                Dict{String,Any}(
                    "id" => "PV1",
                    "componentType" => "PV",
                    "configuration" => Dict{String,Any}("ports" => ["port1"])
                ),
                Dict{String,Any}(
                    "id" => "Motor1", 
                    "componentType" => "motor",
                    "configuration" => Dict{String,Any}("ports" => ["port1"])
                ),
                Dict{String,Any}(
                    "id" => "Bus1",
                    "componentType" => "bus",
                    "configuration" => Dict{String,Any}("ports" => ["port1", "port2", "port3"])
                )
            ],
            "connections" => [
                Dict{String,Any}(
                    "id" => "Cable1",
                    "from" => Dict{String,Any}("componentInstanceId" => "PV1"),
                    "to" => Dict{String,Any}("componentInstanceId" => "Bus1")
                ),
                Dict{String,Any}(
                    "id" => "Cable2",
                    "from" => Dict{String,Any}("componentInstanceId" => "Bus1"),
                    "to" => Dict{String,Any}("componentInstanceId" => "Motor1")
                )
            ]
        );
        
        # Call the dcide-specific solve! function explicitly
        @test _MSS.solve!(test_json) === nothing;
        
        # The function should add STD information to the JSON
        @test haskey(test_json["componentInstances"][1], "std") ||
              haskey(test_json["connections"][1], "std")
    end
end
