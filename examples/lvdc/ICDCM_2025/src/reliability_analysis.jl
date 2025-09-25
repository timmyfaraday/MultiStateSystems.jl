################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Reliability Analysis Module

This module provides functions for analyzing system reliability metrics
including failure rates, availability calculations, and combined failure 
probability analysis for LVDC networks.

## Main Functions
- `calculate_system_failure_rates`: Compute system-level failure rates
- `calculate_load_failure_rates`: Compute load-specific failure rates
- `calculate_combined_failure_probability`: Analyze combined failure scenarios
- `extract_reliability_metrics`: Extract key reliability metrics
"""

using MultiStateSystems
using Unitful

"""
    calculate_system_failure_rates(component_stds, source_stds)

Calculate system-level failure rates from component STDs.

# Arguments
- `component_stds`: Dictionary of component STDs organized by zone/device/feeder
- `source_stds`: Dictionary of source STDs

# Returns
- `Dict`: System failure rates organized by zone and device type
"""
function calculate_system_failure_rates(component_stds, source_stds)
    failure_rates = Dict()
    
    for (zone_key, device_stds) in component_stds
        failure_rates[zone_key] = Dict()
        
        for (device_type, feeder_stds) in device_stds
            # Sum failure rates from all feeders
            total_failure_rate = sum(
                feeder.sprops[3][:h] .+ feeder.sprops[5][:h] 
                for feeder in values(feeder_stds)
            )
            
            # Add source failure rates
            for source in values(source_stds[zone_key]["SSCB"])
                total_failure_rate .+= source.sprops[3][:h] .+ source.sprops[5][:h]
            end
            
            failure_rates[zone_key][device_type] = total_failure_rate
        end
    end
    
    return failure_rates
end

"""
    calculate_load_failure_rates(component_stds, system_failure_rates)

Calculate failure rates specific to each load feeder.

# Arguments
- `component_stds`: Dictionary of component STDs
- `system_failure_rates`: System-level failure rates

# Returns
- `Dict`: Load-specific failure rates
"""
function calculate_load_failure_rates(component_stds, system_failure_rates)
    load_failure_rates = Dict()
    
    for (zone_key, device_stds) in component_stds
        load_failure_rates[zone_key] = Dict()
        
        for (device_type, feeder_stds) in device_stds
            load_failure_rates[zone_key][device_type] = Dict()
            
            for (feeder_id, feeder_std) in feeder_stds
                # Combine system failure rate with feeder-specific failure rate
                feeder_failure_rate = feeder_std.sprops[2][:h]  # Voltage fault state
                system_rate = system_failure_rates[zone_key][device_type]
                
                load_failure_rates[zone_key][device_type][feeder_id] = system_rate .+ feeder_failure_rate
            end
        end
    end
    
    return load_failure_rates
end

"""
    calculate_combined_failure_probability(load_failure_rates, network, zone_key, device_type, feeder1, feeder2, user1, user2, dt)

Calculate combined failure probability for two specific feeders and users.

# Arguments
- `load_failure_rates`: Load-specific failure rates
- `network`: Network object containing user information
- `zone_key`: Zone identifier
- `device_type`: Protection device type
- `feeder1`: First feeder identifier
- `feeder2`: Second feeder identifier  
- `user1`: First user identifier
- `user2`: Second user identifier
- `dt`: Time step for integration

# Returns
- `Vector`: Combined failure probability over time
"""
function calculate_combined_failure_probability(load_failure_rates, network, zone_key, device_type, 
                                              feeder1, feeder2, user1, user2, dt)
    # Extract failure rates
    fail1 = load_failure_rates[zone_key][device_type][feeder1]
    fail2 = load_failure_rates[zone_key][device_type][feeder2]
    
    # Extract user probabilities
    U1 = network[zone_key][device_type].usr[user1][:ugf].prb[1]
    U2 = network[zone_key][device_type].usr[user2][:ugf].prb[1]
    
    # Calculate convolution integrals
    int1 = zeros(length(fail1))*unit(fail1[1])
    int2 = zeros(length(fail2))*unit(fail1[2])
    
    for i in 1:length(fail1)
        for j in 1:i
            int1[i] += fail1[j] * dt * U2[length(fail1)-j+1] * _MSS.weights(i, j)
        end
    end
    
    for i in 1:length(fail2)
        for j in 1:i
            int2[i] += fail2[j] * dt * U1[length(fail2)-j+1] * _MSS.weights(i, j)
        end
    end
    
    return int1 .+ int2
end

"""
    calculate_double_bus_combined_failure(zone1_failure_rates, zone2_failure_rates, network, zone_key, device_type, feeder1, feeder2, user1, user2, dt)

Calculate combined failure probability for double bus configuration.

# Arguments
- `zone1_failure_rates`: Zone 1 failure rates
- `zone2_failure_rates`: Zone 2 failure rates
- `network`: Network object
- `zone_key`: Zone identifier
- `device_type`: Device type
- `feeder1`: Feeder in zone 1
- `feeder2`: Feeder in zone 2
- `user1`: User in zone 1
- `user2`: User in zone 2
- `dt`: Time step

# Returns
- `Vector`: Combined failure probability for double bus system
"""
function calculate_double_bus_combined_failure(zone1_failure_rates, zone2_failure_rates, network, 
                                             zone_key, device_type, feeder1, feeder2, user1, user2, dt)
    # Extract failure rates from different zones
    fail1 = zone1_failure_rates[zone_key][device_type][feeder1]
    fail2 = zone2_failure_rates[zone_key][device_type][feeder2]
    
    # Extract user probabilities
    U1 = network[zone_key][device_type].usr[user1][:ugf].prb[1]
    U2 = network[zone_key][device_type].usr[user2][:ugf].prb[1]
    
    # Calculate cross-zone failure interactions
    int1 = zeros(length(fail1))*unit(fail1[1])
    int2 = zeros(length(fail2))*unit(fail1[2])
    
    for i in 1:length(fail1)
        for j in 1:i
            int1[i] += fail1[j] * dt * U2[length(fail1)-j+1] * _MSS.weights(i, j)
        end
    end
    
    for i in 1:length(fail2)
        for j in 1:i
            int2[i] += fail2[j] * dt * U1[length(fail2)-j+1] * _MSS.weights(i, j)
        end
    end
    
    return int1 .+ int2
end

"""
    analyze_protection_device_performance(fault_probabilities, device_types)

Analyze and compare protection device performance.

# Arguments
- `fault_probabilities`: Fault clearing probabilities for different devices
- `device_types`: List of device types to analyze

# Returns
- `Dict`: Performance comparison metrics
"""
function analyze_protection_device_performance(fault_probabilities, device_types)
    performance = Dict()
    
    for device_type in device_types
        if haskey(fault_probabilities["load"], device_type)
            device_probs = fault_probabilities["load"][device_type]
            
            # Calculate average performance metrics
            avg_clearing_prob = mean(values(device_probs))
            min_clearing_prob = minimum(values(device_probs))
            max_clearing_prob = maximum(values(device_probs))
            
            performance[device_type] = Dict(
                "average_clearing_probability" => avg_clearing_prob,
                "minimum_clearing_probability" => min_clearing_prob,
                "maximum_clearing_probability" => max_clearing_prob,
                "performance_range" => max_clearing_prob - min_clearing_prob
            )
        end
    end
    
    return performance
end

"""
    generate_reliability_report(metrics, performance_analysis)

Generate a comprehensive reliability report.

# Arguments
- `metrics`: Reliability metrics from extract_reliability_metrics
- `performance_analysis`: Device performance analysis

# Returns
- `Dict`: Comprehensive reliability report
"""
function generate_reliability_report(metrics, performance_analysis)
    report = Dict(
        "system_metrics" => metrics,
        "device_performance" => performance_analysis,
        "summary" => Dict(),
        "recommendations" => []
    )
    
    # Generate summary statistics
    all_availabilities = []
    for (zone, devices) in metrics
        for (device, device_metrics) in devices
            push!(all_availabilities, device_metrics["average_availability"])
        end
    end
    
    report["summary"]["overall_average_availability"] = mean(all_availabilities)
    report["summary"]["system_availability_range"] = maximum(all_availabilities) - minimum(all_availabilities)
    
    # Generate recommendations based on performance
    if haskey(performance_analysis, "Fuse")
        fuse_perf = performance_analysis["Fuse"]["average_clearing_probability"]
        if fuse_perf < 0.95
            push!(report["recommendations"], "Consider upgrading fuse protection due to low clearing probability")
        end
    end
    
    if report["summary"]["system_availability_range"] > 0.01
        push!(report["recommendations"], "Investigate availability variance between protection devices")
    end
    
    return report
end
