################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Utility Functions Module

This module contains utility functions and helpers used throughout the LVDC
short-circuit analysis examples. It provides common operations, data processing
functions, and convenience methods.

## Main Functions
- `save_results`: Save analysis results to files
- `load_results`: Load previously saved results
- `print_summary`: Print analysis summary
- `validate_configuration`: Validate system configuration
"""

using MultiStateSystems
using Unitful
using JLD2
using Statistics
using Dates

const _MSS = MultiStateSystems

"""
    save_results(results, filename; format="jld2", filter_keys=nothing)

Save analysis results to a file with optional filtering.

# Arguments
- `results`: Results dictionary to save
- `filename`: Output filename
- `format`: File format ("jld2" or "json", default: "jld2")
- `filter_keys`: Vector of keys to save (default: nothing saves all), or predefined filter strings:
  - "stds_only": Save only solved state transition diagrams
  - "minimal": Save only essential results (STDs + basic metrics)
  - "networks_only": Save only network solutions
  
# Examples
```julia
# Save only solved STDs
save_results(results, "stds.jld2", filter_keys="stds_only")

# Save specific keys
save_results(results, "custom.jld2", filter_keys=["solved_stds", "reliability_metrics"])

# Save all results (default behavior)
save_results(results, "all.jld2")
```
"""
function save_results(results, filename; format="jld2", filter_keys=nothing)
    # Ensure directory exists
    mkpath(dirname(filename))
    
    # Apply filtering if specified
    filtered_results = apply_result_filter(results, filter_keys)
    
    if format == "jld2"
        @save filename filtered_results
        if filter_keys !== nothing
            if isa(filter_keys, String)
                println("Results saved to: $filename (filter: $filter_keys)")
            else
                println("Results saved to: $filename (filtered keys: $(join(filter_keys, ", ")))")
            end
        else
            println("Results saved to: $filename (all data)")
        end
    elseif format == "json"
        # For JSON, we'd need to convert Unitful quantities to dictionaries
        # This is a simplified implementation
        println("JSON export not fully implemented - use JLD2 format")
    else
        error("Unsupported format: $format. Use 'jld2' or 'json'")
    end
end

"""
    apply_result_filter(results, filter_keys)

Apply filtering to results dictionary based on specified keys or predefined filters.

# Arguments
- `results`: Original results dictionary
- `filter_keys`: Filter specification (String for predefined filters, Vector{String} for custom keys, or nothing for no filtering)

# Returns
- `Dict`: Filtered results dictionary
"""
function apply_result_filter(results, filter_keys)
    if filter_keys === nothing
        # No filtering - return all results
        return results
    end
    
    if isa(filter_keys, String)
        # Predefined filter strings
        if filter_keys == "stds_only"
            # Save only solved state transition diagrams
            filtered = Dict()
            if haskey(results, "solved_stds")
                filtered["solved_stds"] = results["solved_stds"]
            end
            # Add metadata for context
            filtered["filter_applied"] = "stds_only"
            filtered["creation_time"] = string(Dates.now())
            return filtered
            
        elseif filter_keys == "minimal"
            # Save STDs + basic reliability metrics
            filtered = Dict()
            minimal_keys = ["solved_stds", "reliability_metrics", "fault_probabilities"]
            for key in minimal_keys
                if haskey(results, key)
                    filtered[key] = results[key]
                end
            end
            filtered["filter_applied"] = "minimal"
            filtered["creation_time"] = string(Dates.now())
            return filtered
            
        elseif filter_keys == "networks_only"
            # Save only network solutions
            filtered = Dict()
            if haskey(results, "networks")
                filtered["networks"] = results["networks"]
            end
            filtered["filter_applied"] = "networks_only"
            filtered["creation_time"] = string(Dates.now())
            return filtered
            
        else
            error("Unknown predefined filter: $filter_keys. Use 'stds_only', 'minimal', or 'networks_only'")
        end
        
    elseif isa(filter_keys, Vector)
        # Custom key filtering
        filtered = Dict()
        for key in filter_keys
            if haskey(results, key)
                filtered[key] = results[key]
            else
                @warn "Key '$key' not found in results"
            end
        end
        filtered["filter_applied"] = "custom"
        filtered["filtered_keys"] = filter_keys
        filtered["creation_time"] = string(Dates.now())
        return filtered
        
    else
        error("filter_keys must be a String (predefined filter), Vector{String} (custom keys), or nothing")
    end
end

"""
    load_results(filename)

Load previously saved analysis results.

# Arguments
- `filename`: Input filename

# Returns
- `Dict`: Loaded results dictionary
"""
function load_results(filename)
    if !isfile(filename)
        error("File not found: $filename")
    end
    
    @load filename results
    println("Results loaded from: $filename")
    return results
end

"""
    print_summary(results; detailed=false)

Print a summary of analysis results.

# Arguments
- `results`: Results dictionary from analysis
- `detailed`: Whether to print detailed information (default: false)
"""
function print_summary(results; detailed=false)
    println("="^60)
    println("LVDC SHORT-CIRCUIT ANALYSIS SUMMARY")
    println("="^60)
    
    # Print configuration information
    if haskey(results, "configuration")
        config = results["configuration"]
        println("Configuration:")
        for (key, value) in config
            println("  $key: $value")
        end
        println()
    end
    
    # Print fault probabilities summary
    if haskey(results, "fault_probabilities")
        println("Fault Clearing Probabilities:")
        fault_probs = results["fault_probabilities"]
        
        for (component_type, devices) in fault_probs
            println("  $component_type:")
            for (device_type, feeders) in devices
                avg_prob = mean(values(feeders))
                min_prob = minimum(values(feeders))
                max_prob = maximum(values(feeders))
                println("    $device_type: avg=$avg_prob, range=[$min_prob, $max_prob]")
            end
        end
        println()
    end
    
    # Print reliability metrics
    if haskey(results, "reliability_metrics") && haskey(results["reliability_metrics"], "user_availability")
        println("User Availability Summary:")
        user_avail = results["reliability_metrics"]["user_availability"]
        
        for (zone, devices) in user_avail
            println("  Zone $zone:")
            for (device, metrics) in devices
                avg_avail = metrics["average_availability"]
                println("    $device: $(round(avg_avail * 100, digits=4))%")
            end
        end
        println()
    end
    
    # Print bridge comparison if available
    if haskey(results, "bridge_comparison")
        bridge_comp = results["bridge_comparison"]
        println("Bridge Configuration Analysis:")
        if haskey(bridge_comp, "summary")
            summary = bridge_comp["summary"]
            improvement = summary["average_availability_improvement"]
            risk_reduction = summary["average_risk_reduction"]
            recommended = summary["bridge_recommended"]
            
            println("  Average availability improvement: $(round(improvement * 100, digits=6))%")
            println("  Average risk reduction: $(round(risk_reduction * 100, digits=6))%")
            println("  Bridge recommended: $recommended")
        end
        println()
    end
    
    if detailed
        print_detailed_summary(results)
    end
    
    println("="^60)
end

"""
    print_detailed_summary(results)

Print detailed analysis results.
"""
function print_detailed_summary(results)
    println("DETAILED ANALYSIS RESULTS:")
    println("-"^40)
    
    # Print detailed failure rates if available
    if haskey(results, "reliability_metrics")
        metrics = results["reliability_metrics"]
        
        for (key, value) in metrics
            if occursin("failure_rates", key)
                println("$key:")
                if isa(value, Dict)
                    for (zone, devices) in value
                        println("  Zone $zone:")
                        for (device, rates) in devices
                            final_rate = rates[end]
                            println("    $device final rate: $final_rate")
                        end
                    end
                end
                println()
            end
        end
    end
    
    # Print combined failure scenarios
    if haskey(results, "combined_failures") || haskey(results, "cross_zone_failures")
        failures_key = haskey(results, "combined_failures") ? "combined_failures" : "cross_zone_failures"
        println("Combined Failure Scenarios:")
        
        for (scenario, failure_prob) in results[failures_key]
            if isa(failure_prob, Vector) && !isempty(failure_prob)
                final_prob = failure_prob[end]
                max_prob = maximum(failure_prob)
                println("  $scenario: final=$final_prob, max=$max_prob")
            end
        end
        println()
    end
end

"""
    validate_configuration(sys_params, int_params)

Validate system and integration parameters.

# Arguments
- `sys_params`: System parameters structure
- `int_params`: Integration parameters structure

# Returns
- `Bool`: True if configuration is valid

# Throws
- `ArgumentError`: If configuration is invalid
"""
function validate_configuration(sys_params, int_params)
    # Validate system parameters
    if sys_params.V_DC <= 0
        throw(ArgumentError("DC voltage must be positive"))
    end
    
    if sys_params.I_max <= 0
        throw(ArgumentError("Maximum current must be positive"))
    end
    
    if sys_params.V_min >= sys_params.V_DC
        throw(ArgumentError("Minimum voltage must be less than DC voltage"))
    end
    
    if sys_params.C_b <= 0
        throw(ArgumentError("Bus capacitance must be positive"))
    end
    
    if any(p -> p < 0 || p > 1, sys_params.P_zone)
        throw(ArgumentError("Zone probabilities must be between 0 and 1"))
    end
    
    # Validate integration parameters
    if ustrip(int_params.t_max) <= 0
        throw(ArgumentError("Maximum simulation time must be positive"))
    end
    
    if int_params.n <= 0
        throw(ArgumentError("Number of integration points must be positive"))
    end
    
    if ustrip(int_params.tsim) <= 0
        throw(ArgumentError("Simulation time must be positive"))
    end
    
    if ustrip(int_params.dt) <= 0
        throw(ArgumentError("Time step must be positive"))
    end
    
    println("Configuration validation passed ✓")
    return true
end

"""
    calculate_performance_metrics(availability_data)

Calculate various performance metrics from availability data.

# Arguments
- `availability_data`: Dictionary of availability data

# Returns
- `Dict`: Performance metrics including MTBF, downtime, etc.
"""
function calculate_performance_metrics(availability_data)
    metrics = Dict()
    
    for (zone, devices) in availability_data
        metrics[zone] = Dict()
        
        for (device, device_metrics) in devices
            if haskey(device_metrics, "average_availability")
                availability = device_metrics["average_availability"]
                
                # Calculate metrics
                downtime_fraction = 1 - availability
                annual_downtime_hours = downtime_fraction * 8760  # hours per year
                annual_downtime_minutes = annual_downtime_hours * 60
                
                # Simplified MTBF calculation (assuming constant failure rate)
                mtbf_hours = availability > 0 ? 8760 / downtime_fraction : Inf
                
                metrics[zone][device] = Dict(
                    "availability_percent" => availability * 100,
                    "downtime_fraction" => downtime_fraction,
                    "annual_downtime_hours" => annual_downtime_hours,
                    "annual_downtime_minutes" => annual_downtime_minutes,
                    "mtbf_hours" => mtbf_hours,
                    "reliability_class" => classify_reliability(availability)
                )
            end
        end
    end
    
    return metrics
end

"""
    classify_reliability(availability)

Classify reliability level based on availability.

# Arguments
- `availability`: Availability fraction (0-1)

# Returns
- `String`: Reliability classification
"""
function classify_reliability(availability)
    if availability >= 0.9999
        return "Excellent (>99.99%)"
    elseif availability >= 0.999
        return "Very Good (>99.9%)"
    elseif availability >= 0.99
        return "Good (>99%)"
    elseif availability >= 0.95
        return "Acceptable (>95%)"
    else
        return "Poor (<95%)"
    end
end

"""
    create_results_directory(base_path="results")

Create a timestamped results directory.

# Arguments
- `base_path`: Base directory for results (default: "results")

# Returns
- `String`: Path to created directory
"""
function create_results_directory(base_path="results")
    timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
    results_dir = joinpath(base_path, "analysis_$timestamp")
    mkpath(results_dir)
    return results_dir
end

"""
    export_to_csv(data, filename)

Export numerical data to CSV format.

# Arguments
- `data`: Dictionary or array of numerical data
- `filename`: Output CSV filename
"""
function export_to_csv(data, filename)
    # This would require the CSV.jl package for full implementation
    # For now, provide a basic implementation for simple data
    
    mkpath(dirname(filename))
    
    if isa(data, Dict)
        # Simple CSV export for dictionary data
        open(filename, "w") do io
            # Write header
            keys_list = collect(keys(data))
            println(io, join(keys_list, ","))
            
            # Write data (assuming all values are scalars or same-length vectors)
            first_value = first(values(data))
            if isa(first_value, Number)
                # Single row of data
                values_list = [data[k] for k in keys_list]
                println(io, join(values_list, ","))
            elseif isa(first_value, AbstractVector)
                # Multiple rows
                n_rows = length(first_value)
                for i in 1:n_rows
                    values_list = [data[k][i] for k in keys_list]
                    println(io, join(values_list, ","))
                end
            end
        end
        
        println("Data exported to: $filename")
    else
        println("CSV export only supports Dictionary data currently")
    end
end

"""
    save_solved_stds(solved_stds, filename)

Save only the solved state transition diagrams to a JLD2 file.

# Arguments
- `solved_stds`: Dictionary of solved state transition diagrams
- `filename`: Output filename

# Examples
```julia
save_solved_stds(results["solved_stds"], "results/single_bus_stds.jld2")
```
"""
function save_solved_stds(solved_stds, filename)
    # Ensure directory exists
    mkpath(dirname(filename))
    
    # Create a minimal results structure with only STDs
    std_results = Dict(
        "solved_stds" => solved_stds,
        "creation_time" => string(Dates.now()),
        "data_type" => "solved_stds_only"
    )
    
    @save filename std_results
    println("Solved STDs saved to: $filename")
end

"""
    load_solved_stds(filename)

Load previously solved state transition diagrams from a JLD2 file.

# Arguments
- `filename`: Input filename

# Returns
- `Dict`: Dictionary of solved state transition diagrams

# Examples
```julia
solved_stds = load_solved_stds("results/single_bus_stds.jld2")
```
"""
function load_solved_stds(filename)
    if !isfile(filename)
        error("STD file not found: $filename")
    end
    
    @load filename std_results
    
    if !haskey(std_results, "solved_stds")
        error("File does not contain solved STDs: $filename")
    end
    
    println("Solved STDs loaded from: $filename")
    return std_results["solved_stds"]
end
