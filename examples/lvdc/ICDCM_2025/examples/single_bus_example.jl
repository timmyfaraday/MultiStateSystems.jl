################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Single Bus Example

This module implements the single bus LVDC network analysis example,
demonstrating short-circuit analysis and reliability calculations for
a radial LVDC distribution system.

## Functions
- `single_bus_analysis`: Main analysis function
- `compute_fault_probabilities`: Calculate fault clearing probabilities  
- `build_single_bus_networks`: Construct network topologies
- `analyze_single_bus_reliability`: Perform reliability analysis
"""

using MultiStateSystems
using Unitful
using JLD2

# Load required modules
include("../src/fault_analysis.jl")
include("../src/std_generation.jl")
include("../src/network_builder.jl")
include("../src/reliability_analysis.jl")

"""
    single_bus_analysis(sys_params, int_params; process_type="SemiMarkov", load_stds_from=nothing)

Perform complete single bus LVDC network analysis.

# Arguments
- `sys_params`: System parameters structure
- `int_params`: Integration parameters structure
- `process_type`: "SemiMarkov" or "Markov" (default: "SemiMarkov")
- `load_stds_from`: Path to pre-solved STDs file (default: nothing, solves from scratch)

# Returns
- `Dict`: Complete analysis results including networks, STDs, and metrics
"""
function single_bus_analysis(sys_params, int_params; process_type="SemiMarkov", load_stds_from=nothing)
    println("Starting single bus analysis...")
    
    results = Dict()
    
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
        
        # Step 1a: Calculate fault clearing probabilities
        println("    Computing fault clearing probabilities...")
        fault_probabilities = compute_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)
        converter_probabilities = compute_converter_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)
        
        results["fault_probabilities"] = fault_probabilities
        results["converter_probabilities"] = converter_probabilities
        
        # Step 1b: Generate and solve STDs for all protection zones
        println("    Generating state transition diagrams...")
        solved_stds = Dict()
        
        for P_zone in sys_params.P_zone
            println("      Processing zone with P_zone = $P_zone")
            
            # Generate STDs
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
            
            # Solve STDs
            process_class = process_type == "SemiMarkov" ? SemiMarkovProcess() : MarkovProcess()
            solve_component_stds!(load_stds, process_class, int_params.tsim, int_params.dt)
            solve_component_stds!(source_stds, process_class, int_params.tsim, int_params.dt)
            
            # Create reduced representations if using SemiMarkov
            if process_type == "SemiMarkov"
                load_stds = create_reduced_stds(Dict("load" => load_stds))["load"]
                source_stds = create_reduced_stds(Dict("source" => source_stds))["source"]
            end
            
            solved_stds[P_zone] = Dict("loads" => load_stds, "sources" => source_stds)
        end
        
        results["solved_stds"] = solved_stds
        println("  ✓ STDs solved successfully")
    end
    
    # Step 2: Build and solve networks (always done regardless of STD source)
    println("  Building network topologies...")
    networks = build_single_bus_networks(results["solved_stds"], get_ac_grid_availability())
    
    println("  Solving networks...")
    solve_multiple_networks!(networks)
    
    results["networks"] = networks
    
    # Step 3: Reliability analysis
    println("  Performing reliability analysis...")
    reliability_metrics = analyze_single_bus_reliability(results["solved_stds"], networks, int_params.tsim)
    results["reliability_metrics"] = reliability_metrics
    
    # Step 4: Calculate combined failure probabilities (example)
    println("  Computing combined failure scenarios...")
    combined_failures = calculate_single_bus_combined_failures(results["solved_stds"], networks, int_params.dt)
    results["combined_failures"] = combined_failures
    
    println("Single bus analysis completed successfully!")
    return results
end

"""
    compute_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)

Compute fault clearing probabilities for all devices and cable configurations.

# Arguments
- `sys_params`: System parameters
- `int_params`: Integration parameters  
- `device_cable_mapping`: Device to cable mapping
- `protection_devices`: Protection device specifications

# Returns
- `Dict`: Fault clearing probabilities organized by component type, device, and feeder
"""
function compute_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)
    fault_probs = Dict()
    
    for (component_type, device_configs) in device_cable_mapping
        fault_probs[component_type] = Dict()
        
        for (device_type, cable_configs) in device_configs
            fault_probs[component_type][device_type] = Dict()
            
            # Handle combined device types (e.g., "Fuse_MCCB")
            if occursin("_", device_type)
                device_a, device_b = split(device_type, "_")
                
                for (cable_id, cable_lengths) in cable_configs
                    # Calculate probability vectors for both devices
                    # This is crucial: we need vectors to compare probabilities at each cable location
                    prob_vector_a = calculate_device_fault_probability_vector(device_a, cable_lengths, sys_params, int_params, protection_devices)
                    prob_vector_b = calculate_device_fault_probability_vector(device_b, cable_lengths, sys_params, int_params, protection_devices)
                    
                    # Take maximum at each location (best protection), then mean over all locations
                    fault_probs[component_type][device_type][cable_id] = mean(max.(prob_vector_a, prob_vector_b))
                end
            else
                for (cable_id, cable_lengths) in cable_configs
                    fault_probs[component_type][device_type][cable_id] = 
                        calculate_device_fault_probability(device_type, cable_lengths, sys_params, int_params, protection_devices)
                end
            end
        end
    end
    
    return fault_probs
end

"""
    calculate_device_fault_probability(device_type, cable_lengths, sys_params, int_params, protection_devices)

Calculate fault clearing probability for a specific device and cable configuration.
"""
function calculate_device_fault_probability(device_type, cable_lengths, sys_params, int_params, protection_devices)
    device = protection_devices[device_type]
    
    if device.type == "CB"
        return calculate_fault_clearing_probability(
            cable_lengths, sys_params, int_params, (device.μ, device.λ); return_mean=true
        )
    elseif device.type == "Fuse"
        return calculate_fuse_clearing_probability(
            cable_lengths, sys_params, int_params, device.λ; return_mean=true
        )
    else
        error("Unknown device type: $(device.type)")
    end
end

"""
    calculate_device_fault_probability_vector(device_type, cable_lengths, sys_params, int_params, protection_devices)

Calculate fault clearing probability vector for a specific device and cable configuration.
Returns the full probability vector at each cable location instead of the mean.

This function is essential for combined protection devices where probabilities must be 
compared at each cable location before taking the mean. For example, with "Fuse_MCCB" 
protection, we need:
1. Calculate probability vector for Fuse at each cable location
2. Calculate probability vector for MCCB at each cable location  
3. Take maximum probability at each location (best protection)
4. Calculate mean of the resulting maximum probability vector

# Arguments
- `device_type`: Protection device type ("SSCB", "MCCB", "Fuse")
- `cable_lengths`: Cable length range or array
- `sys_params`: System parameters structure
- `int_params`: Integration parameters structure  
- `protection_devices`: Protection device specifications

# Returns
- `Vector` or `Number`: Probability vector for range inputs, single value for scalar inputs
"""
function calculate_device_fault_probability_vector(device_type, cable_lengths, sys_params, int_params, protection_devices)
    device = protection_devices[device_type]
    
    if device.type == "CB"
        return calculate_fault_clearing_probability(
            cable_lengths, sys_params, int_params, (device.μ, device.λ); return_mean=false
        )
    elseif device.type == "Fuse"
        return calculate_fuse_clearing_probability(
            cable_lengths, sys_params, int_params, device.λ; return_mean=false
        )
    else
        error("Unknown device type: $(device.type)")
    end
end

"""
    compute_converter_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)

Compute worst-case fault clearing probabilities for converter protection.
"""
function compute_converter_fault_probabilities(sys_params, int_params, device_cable_mapping, protection_devices)
    converter_probs = Dict()
    
    for (component_type, device_configs) in device_cable_mapping
        converter_probs[component_type] = Dict()
        
        for (device_type, cable_configs) in device_configs
            converter_probs[component_type][device_type] = Dict()
            
            if occursin("_", device_type)
                device_a, device_b = split(device_type, "_")
                
                for (cable_id, cable_lengths) in cable_configs
                    # Use maximum cable length for worst case
                    max_length = maximum(cable_lengths)
                    prob_vector_a = calculate_device_fault_probability_vector(device_a, max_length, sys_params, int_params, protection_devices)
                    prob_vector_b = calculate_device_fault_probability_vector(device_b, max_length, sys_params, int_params, protection_devices)
                    
                    # Handle single length case (returns single value, not vector)
                    if isa(prob_vector_a, Number) && isa(prob_vector_b, Number)
                        converter_probs[component_type][device_type][cable_id] = max(prob_vector_a, prob_vector_b)
                    else
                        converter_probs[component_type][device_type][cable_id] = mean(max.(prob_vector_a, prob_vector_b))
                    end
                end
            else
                for (cable_id, cable_lengths) in cable_configs
                    max_length = maximum(cable_lengths)
                    prob_vector = calculate_device_fault_probability_vector(device_type, max_length, sys_params, int_params, protection_devices)
                    
                    # Handle single length case (returns single value, not vector)
                    if isa(prob_vector, Number)
                        converter_probs[component_type][device_type][cable_id] = prob_vector
                    else
                        converter_probs[component_type][device_type][cable_id] = mean(prob_vector)
                    end
                end
            end
        end
    end
    
    return converter_probs
end

"""
    build_single_bus_networks(solved_stds, ac_availability)

Build single bus networks for all protection zones and device types.
"""
function build_single_bus_networks(solved_stds, ac_availability)
    networks = Dict()
    
    for (P_zone, stds) in solved_stds
        networks[P_zone] = Dict()
        
        for (device_type, load_feeders) in stds["loads"]
            # Get corresponding source STDs
            source_feeders = stds["sources"]["SSCB"]  # Assuming SSCB for sources
            
            # Create AC grid STD
            time_vector = first(values(load_feeders)).props[:time]
            ac_grid_std = create_ac_grid_std(time_vector, ac_availability)
            
            # Build network
            networks[P_zone][device_type] = build_single_bus_network(
                load_feeders, source_feeders, ac_grid_std
            )
        end
    end
    
    return networks
end

"""
    analyze_single_bus_reliability(solved_stds, networks, time_span)

Perform comprehensive reliability analysis for single bus networks.
"""
function analyze_single_bus_reliability(solved_stds, networks, time_span)
    metrics = Dict()

    # Calculate system failure rates
    for (P_zone, stds) in solved_stds
        metrics["failure_rates_$(P_zone)"] = calculate_system_failure_rates(
            Dict(P_zone => stds["loads"]), 
            Dict(P_zone => stds["sources"])
        )
    end
    
    return metrics
end

"""
    calculate_single_bus_combined_failures(solved_stds, networks, dt)

Calculate combined failure probabilities for example feeder pairs.
"""
function calculate_single_bus_combined_failures(solved_stds, networks, dt)
    combined_failures = Dict()
    
    # Example: Calculate combined failure for C2 and C4 feeders with MCCB protection
    P_zone = 1.0  # Use specific zone
    device_type = "MCCB"
    
    if haskey(solved_stds, P_zone) && haskey(solved_stds[P_zone]["loads"], device_type)
        # Calculate system failure rates
        system_rates = calculate_system_failure_rates(
            Dict(P_zone => solved_stds[P_zone]["loads"]),
            Dict(P_zone => solved_stds[P_zone]["sources"])
        )
        
        # Calculate load-specific failure rates  
        load_rates = calculate_load_failure_rates(
            Dict(P_zone => solved_stds[P_zone]["loads"]),
            system_rates
        )
        
        # Calculate combined failure for specific feeders
        if haskey(load_rates[P_zone][device_type], "C2") && haskey(load_rates[P_zone][device_type], "C4")
            combined_failures["C2_C4_MCCB"] = calculate_combined_failure_probability(
                load_rates, networks, P_zone, device_type, "C2", "C4", 3, 5, ustrip(dt)
            )
        end
    end
    
    return combined_failures
end
