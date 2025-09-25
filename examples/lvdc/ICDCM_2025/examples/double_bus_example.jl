################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Double Bus Example

This module implements the double bus LVDC network analysis example,
demonstrating short-circuit analysis and reliability calculations for
a double bus LVDC distribution system with optional bridge connections.

## Functions
- `double_bus_analysis`: Main analysis function
- `build_double_bus_networks`: Construct double bus network topologies
- `analyze_double_bus_reliability`: Perform reliability analysis
- `calculate_double_bus_combined_failures`: Analyze cross-zone failure scenarios
"""

using MultiStateSystems
using Unitful
using JLD2

# Load required modules
include("../src/fault_analysis.jl")
include("../src/std_generation.jl")
include("../src/network_builder.jl")
include("../src/reliability_analysis.jl")
include("single_bus_example.jl")  # Reuse single bus functions

"""
    double_bus_analysis(sys_params, int_params; process_type="SemiMarkov", include_bridge=true, load_stds_from=nothing)

Perform complete double bus LVDC network analysis.

# Arguments
- `sys_params`: System parameters structure
- `int_params`: Integration parameters structure
- `process_type`: "SemiMarkov" or "Markov" (default: "SemiMarkov")
- `include_bridge`: Whether to include bridge cable between buses (default: true)
- `load_stds_from`: Path to pre-solved STDs file (default: nothing, solves from scratch)

# Returns
- `Dict`: Complete analysis results including networks, STDs, and metrics
"""
function double_bus_analysis(sys_params, int_params; process_type="SemiMarkov", include_bridge=true, load_stds_from=nothing)
    println("Starting double bus analysis...")
    
    results = Dict()
    results["configuration"] = Dict("include_bridge" => include_bridge, "process_type" => process_type)
    
    # Step 1: Get or load solved STDs
    if load_stds_from !== nothing
        println("  Loading pre-solved STDs from: $load_stds_from")
        solved_stds = load_solved_stds(load_stds_from)
        results["solved_stds"] = solved_stds
        println("  ✓ STDs loaded successfully")
    else
        println("  Solving STDs from scratch...")
        
        # Get configuration data
        device_cable_mapping = get_device_cable_mapping()
        protection_devices = get_protection_devices()
        zone_config = get_network_zones()
        
        # Step 1a: Calculate fault clearing probabilities (reuse from single bus)
        println("    Computing fault clearing probabilities...")
        fault_probabilities = compute_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)
        converter_probabilities = compute_converter_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)
        
        results["fault_probabilities"] = fault_probabilities
        results["converter_probabilities"] = converter_probabilities
        
        # Step 1b: Generate and solve STDs including bridge if requested
        println("    Generating state transition diagrams...")
        solved_stds = Dict()
    
    for P_zone in sys_params.P_zone
        println("    Processing zone with P_zone = $P_zone")
        
        # Generate load and source STDs (same as single bus)
        load_stds = generate_component_stds(
            Dict("load" => device_cable_mapping["load"]),
            Dict("load" => fault_probabilities["load"]),
            sys_params.λᶜ;
            converter_probabilities = Dict("load" => converter_probabilities["load"]),
            process_type = process_type,
            P_zone = P_zone
        )["load"]
        
        source_stds = generate_component_stds(
            Dict("source" => device_cable_mapping["source"]),
            Dict("source" => fault_probabilities["source"]),
            sys_params.λᶜ;
            converter_probabilities = Dict("source" => converter_probabilities["source"]),
            process_type = process_type,
            P_zone = P_zone
        )["source"]
        
        # Generate bridge STDs if requested
        bridge_stds = nothing
        if include_bridge
            bridge_stds = generate_component_stds(
                Dict("bridge" => device_cable_mapping["bridge"]),
                Dict("bridge" => fault_probabilities["bridge"]),
                sys_params.λᶜ;
                process_type = process_type,
                P_zone = P_zone
            )["bridge"]
        end
        
        # Solve all STDs
        process_class = process_type == "SemiMarkov" ? SemiMarkovProcess() : MarkovProcess()
        solve_component_stds!(load_stds, process_class, int_params.tsim, int_params.dt)
        solve_component_stds!(source_stds, process_class, int_params.tsim, int_params.dt)
        
        if bridge_stds !== nothing
            solve_component_stds!(bridge_stds, process_class, int_params.tsim, int_params.dt)
        end
        
        # Create reduced representations if using SemiMarkov
        if process_type == "SemiMarkov"
            load_stds = create_reduced_stds(Dict("load" => load_stds))["load"]
            source_stds = create_reduced_stds(Dict("source" => source_stds))["source"]
            if bridge_stds !== nothing
                bridge_stds = create_reduced_stds(Dict("bridge" => bridge_stds))["bridge"]
            end
        end
        
        solved_stds[P_zone] = Dict(
            "loads" => load_stds, 
            "sources" => source_stds,
            "bridge" => bridge_stds
        )
    end
    
    results["solved_stds"] = solved_stds
        println("  ✓ STDs solved successfully")
    end
    
    # Step 2: Build and solve double bus networks (always done regardless of STD source)
    println("  Building double bus network topologies...")
    zone_config = get_network_zones()  # Get zone configuration for network building
    networks = build_double_bus_networks(results["solved_stds"], zone_config, get_ac_grid_availability(), include_bridge)
    
    println("  Solving networks...")
    solve_multiple_networks!(networks)
    
    results["networks"] = networks
    
    # Step 3: Reliability analysis
    println("  Performing reliability analysis...")
    reliability_metrics = analyze_double_bus_reliability(results["solved_stds"], networks, zone_config, int_params.tsim)
    results["reliability_metrics"] = reliability_metrics
    
    # Step 4: Calculate cross-zone failure scenarios
    println("  Computing cross-zone failure scenarios...")
    cross_zone_failures = calculate_double_bus_combined_failures(results["solved_stds"], networks, zone_config, int_params.dt)
    results["cross_zone_failures"] = cross_zone_failures

    println("Double bus analysis completed successfully!")
    return results
end

"""
    build_double_bus_networks(solved_stds, zone_config, ac_availability, include_bridge)

Build double bus networks for all protection zones and device types.
"""
function build_double_bus_networks(solved_stds, zone_config, ac_availability, include_bridge)
    networks = Dict()
    
    for (P_zone, stds) in solved_stds
        networks[P_zone] = Dict()
        
        for (device_type, load_feeders) in stds["loads"]
            # Get corresponding source and bridge STDs
            source_feeders = stds["sources"]["SSCB"]
            bridge_feeders = include_bridge ? stds["bridge"] : nothing
            
            # Create AC grid STD
            time_vector = first(values(load_feeders)).props[:time]
            ac_grid_std = create_ac_grid_std(time_vector, ac_availability)
            
            # Build double bus network
            networks[P_zone][device_type] = build_double_bus_network(
                load_feeders, source_feeders, ac_grid_std, zone_config; 
                bridge_stds = bridge_feeders
            )
        end
    end
    
    return networks
end

"""
    analyze_double_bus_reliability(solved_stds, networks, zone_config, time_span)

Perform comprehensive reliability analysis for double bus networks.
"""
function analyze_double_bus_reliability(solved_stds, networks, zone_config, time_span)
    metrics = Dict()

    # Calculate zone-specific failure rates
    for (P_zone, stds) in solved_stds
        # Zone 1 failure rates
        zone_1_loads = Dict(k => v for (k, v) in stds["loads"] if k in ["SSCB", "MCCB", "Fuse"])
        zone_1_filtered = Dict()
        for (device_type, feeders) in zone_1_loads
            zone_1_filtered[device_type] = Dict(k => v for (k, v) in feeders if k in zone_config["Zone_1"])
        end
        
        metrics["zone_1_failure_rates_$(P_zone)"] = calculate_system_failure_rates(
            Dict(P_zone => zone_1_filtered),
            Dict(P_zone => Dict("SSCB" => Dict("S1" => stds["sources"]["SSCB"]["S1"])))
        )
        
        # Zone 2 failure rates
        zone_2_filtered = Dict()
        for (device_type, feeders) in zone_1_loads
            zone_2_filtered[device_type] = Dict(k => v for (k, v) in feeders if k in zone_config["Zone_2"])
        end
        
        metrics["zone_2_failure_rates_$(P_zone)"] = calculate_system_failure_rates(
            Dict(P_zone => zone_2_filtered),
            Dict(P_zone => Dict("SSCB" => Dict("S2" => stds["sources"]["SSCB"]["S2"])))
        )
    end
    
    # Calculate cross-zone dependencies
    metrics["cross_zone_analysis"] = analyze_cross_zone_dependencies(solved_stds, zone_config)
    
    return metrics
end

"""
    analyze_cross_zone_dependencies(solved_stds, zone_config)

Analyze how failures in one zone affect the other zone.
"""
function analyze_cross_zone_dependencies(solved_stds, zone_config)
    dependencies = Dict()
    
    for (P_zone, stds) in solved_stds
        # Analyze how Zone 1 failures affect Zone 2 and vice versa
        zone_1_impact = 0.0
        zone_2_impact = 0.0
        
        for (device_type, feeders) in stds["loads"]
            for (feeder_id, std) in feeders
                # Check for severe unavailability states that propagate
                for state in values(std.sprops)
                    if occursin("U3", state[:name])  # Severe failures
                        final_prob = state[:prob][end]
                        if feeder_id in zone_config["Zone_1"]
                            zone_1_impact += final_prob * P_zone
                        elseif feeder_id in zone_config["Zone_2"]
                            zone_2_impact += final_prob * P_zone
                        end
                    end
                end
            end
        end
        
        dependencies[P_zone] = Dict(
            "zone_1_impact_on_zone_2" => zone_1_impact,
            "zone_2_impact_on_zone_1" => zone_2_impact,
            "total_cross_zone_risk" => zone_1_impact + zone_2_impact
        )
    end
    
    return dependencies
end

"""
    calculate_double_bus_combined_failures(solved_stds, networks, zone_config, dt)

Calculate combined failure probabilities for cross-zone scenarios.
"""
function calculate_double_bus_combined_failures(solved_stds, networks, zone_config, dt)
    combined_failures = Dict()
    
    # Example: Calculate combined failure for C2 (Zone 1) and C4 (Zone 2) with SSCB protection
    P_zone = 1.0
    device_type = "SSCB"
    
    if haskey(solved_stds, P_zone) && haskey(solved_stds[P_zone]["loads"], device_type)
        # Calculate zone-specific failure rates
        zone_1_loads = Dict(k => v for (k, v) in solved_stds[P_zone]["loads"][device_type] if k in zone_config["Zone_1"])
        zone_2_loads = Dict(k => v for (k, v) in solved_stds[P_zone]["loads"][device_type] if k in zone_config["Zone_2"])
        
        # Calculate system failure rates for each zone
        zone_1_system_rates = calculate_system_failure_rates(
            Dict(P_zone => Dict(device_type => zone_1_loads)),
            Dict(P_zone => Dict("SSCB" => Dict("S1" => solved_stds[P_zone]["sources"]["SSCB"]["S1"])))
        )
        
        zone_2_system_rates = calculate_system_failure_rates(
            Dict(P_zone => Dict(device_type => zone_2_loads)),
            Dict(P_zone => Dict("SSCB" => Dict("S2" => solved_stds[P_zone]["sources"]["SSCB"]["S2"])))
        )
        
        # Calculate load-specific failure rates
        zone_1_load_rates = calculate_load_failure_rates(
            Dict(P_zone => Dict(device_type => zone_1_loads)),
            zone_1_system_rates
        )
        
        zone_2_load_rates = calculate_load_failure_rates(
            Dict(P_zone => Dict(device_type => zone_2_loads)),
            zone_2_system_rates
        )
        
        # Calculate cross-zone combined failure
        if haskey(zone_1_load_rates[P_zone][device_type], "C2") && haskey(zone_2_load_rates[P_zone][device_type], "C4")
            combined_failures["C2_C4_cross_zone"] = calculate_double_bus_combined_failure(
                zone_1_load_rates, zone_2_load_rates, networks, P_zone, device_type, "C2", "C4", 4, 6, ustrip(dt)
            )
        end
    end
    
    return combined_failures
end

