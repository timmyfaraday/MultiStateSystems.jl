################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    NetworkAnalysis

Module for analyzing protection device reliability in LVDC networks.
Contains functions for calculating availability, battery backup analysis, and bus reliability.
"""

using MultiStateSystems
using Unitful

const _MSS = MultiStateSystems

"""
    BatteryAnalysis

Structure to hold battery system analysis results and parameters.
"""
struct BatteryAnalysis
    backup_time::typeof(1.0u"hr")
    availability_over_time::Vector{Float64}
    repair_distribution::_MSS.AbstractLogNormal
    
    function BatteryAnalysis(backup_time, availability, repair_mean, repair_std)
        # Use MultiStateSystems' LogNormal distribution with units stripped for compatibility
        repair_dist = _MSS.LogNormal(ustrip(u"d", repair_mean), ustrip(u"d", repair_std))
        new(backup_time, availability, repair_dist)
    end
end

"""
    calculate_battery_availability(
        failure_prob::Float64, 
        backup_time, 
        repair_distribution::_MSS.AbstractLogNormal
    ) -> Float64

Calculate probability that battery backup will suffice until repair is complete.

# Arguments
- `failure_prob::Float64`: Probability of primary power failure
- `backup_time`: Battery backup duration
- `repair_distribution::_MSS.AbstractLogNormal`: Repair time distribution

# Returns
- `Float64`: Probability that battery backup is sufficient
"""
function calculate_battery_availability(
    failure_prob::Float64, 
    backup_time, 
    repair_distribution::_MSS.AbstractLogNormal
)
    # Probability that repair takes longer than battery backup time
    backup_time_days = ustrip(u"d", backup_time)
    prob_repair_exceeds_backup = _MSS.ccdf(repair_distribution, backup_time_days)
    
    # Probability that battery backup is insufficient
    prob_backup_insufficient = failure_prob * prob_repair_exceeds_backup
    
    return 1.0 - prob_backup_insufficient
end

"""
    analyze_battery_system(
        network_availability::Vector{Float64},
        backup_time,
        repair_mean,
        repair_std,
        time_step
    ) -> BatteryAnalysis

Analyze battery backup system performance over time.

# Arguments
- `network_availability::Vector{Float64}`: Network availability over time
- `backup_time`: Battery backup duration
- `repair_mean`: Mean repair time
- `repair_std`: Standard deviation of repair time
- `time_step`: Simulation time step

# Returns
- `BatteryAnalysis`: Complete battery analysis results
"""
function analyze_battery_system(
    network_availability::Vector{Float64},
    backup_time,
    repair_mean,
    repair_std,
    time_step
)
    # Create MultiStateSystems LogNormal distribution for repair times
    repair_dist = _MSS.LogNormal(ustrip(u"d", repair_mean), ustrip(u"d", repair_std))
    backup_steps = Int64(round(backup_time / time_step))
    
    # Initialize availability vector
    battery_availability = Float64[]
    
    # Perfect availability during initial battery backup period
    append!(battery_availability, ones(backup_steps))
    
    # Calculate availability for remaining time
    remaining_steps = length(network_availability) - backup_steps
    for i in 1:remaining_steps
        failure_prob = 1.0 - network_availability[backup_steps + i]
        availability = calculate_battery_availability(failure_prob, backup_time, repair_dist)
        push!(battery_availability, availability)
    end
    
    return BatteryAnalysis(backup_time, battery_availability, repair_dist)
end

"""
    BusAvailabilityResult

Structure to hold bus availability analysis results.
"""
struct BusAvailabilityResult
    protection_device::String
    availability_prob::Vector{Float64}
    unavailability_prob::Vector{Float64}
    time_vector::Vector
    excluded_feeders::Vector{String}
    
    function BusAvailabilityResult(device, avail, unavail, time, excluded=[])
        new(device, avail, unavail, time, excluded)
    end
end

"""
    calculate_bus_availability(
        std_data::Dict,
        time_vector::Vector,
        config;
        excluded_feeders::Vector{String} = String[]
    ) -> Dict{String, BusAvailabilityResult}

Calculate DC bus availability for each protection device type.

# Arguments
- `std_data::Dict`: Processed STD data
- `time_vector::Vector`: Time vector for analysis
- `config`: Configuration object
- `excluded_feeders::Vector{String}`: Feeders to exclude from analysis

# Returns
- `Dict{String, BusAvailabilityResult}`: Bus availability results by protection device
"""
function calculate_bus_availability(
    std_data::Dict,
    time_vector::Vector,
    config;
    excluded_feeders::Vector{String} = String[]
)
    results = Dict{String, BusAvailabilityResult}()
    
    for device_type in config.network.protection_devices
        if !haskey(std_data, device_type)
            @warn "Protection device $device_type not found in STD data"
            continue
        end
        
        # Calculate total failure probability
        total_failure_prob = zeros(length(time_vector))
        
        for (feeder_id, feeder_std) in std_data[device_type]
            # Skip excluded feeders
            if feeder_id in excluded_feeders
                continue
            end
            
            try
                # Get failure probability (assuming state 3 is failure state)
                failure_prob = _MSS.get_sprop(feeder_std, :prob)[3]
                total_failure_prob .+= failure_prob
            catch e
                @warn "Failed to extract failure probability for $device_type/$feeder_id: $e"
                continue
            end
        end
        
        # Calculate availability and unavailability
        unavailability = min.(total_failure_prob, 1.0)  # Cap at 1.0
        availability = 1.0 .- unavailability
        
        results[device_type] = BusAvailabilityResult(
            device_type,
            availability,
            unavailability,
            time_vector,
            excluded_feeders
        )
    end
    
    return results
end

"""
    calculate_selective_bus_availability(
        std_data::Dict,
        time_vector::Vector,
        config,
        included_feeders::Vector{String}
    ) -> Dict{String, BusAvailabilityResult}

Calculate DC bus availability considering only specific feeders.

# Arguments
- `std_data::Dict`: Processed STD data
- `time_vector::Vector`: Time vector for analysis
- `config`: Configuration object
- `included_feeders::Vector{String}`: Feeders to include in analysis

# Returns
- `Dict{String, BusAvailabilityResult}`: Bus availability results by protection device
"""
function calculate_selective_bus_availability(
    std_data::Dict,
    time_vector::Vector,
    config,
    included_feeders::Vector{String}
)
    results = Dict{String, BusAvailabilityResult}()
    
    for device_type in config.network.protection_devices
        if !haskey(std_data, device_type)
            @warn "Protection device $device_type not found in STD data"
            continue
        end
        
        # Calculate total failure probability for included feeders only
        total_failure_prob = zeros(length(time_vector))
        
        for feeder_id in included_feeders
            if !haskey(std_data[device_type], feeder_id)
                @warn "Feeder $feeder_id not found for device $device_type"
                continue
            end
            
            try
                # Get failure probability
                failure_prob = _MSS.get_sprop(std_data[device_type][feeder_id], :prob)[3]
                total_failure_prob .+= failure_prob
            catch e
                @warn "Failed to extract failure probability for $device_type/$feeder_id: $e"
                continue
            end
        end
        
        # Calculate availability and unavailability
        unavailability = min.(total_failure_prob, 1.0)  # Cap at 1.0
        availability = 1.0 .- unavailability
        
        results[device_type] = BusAvailabilityResult(
            device_type,
            availability,
            unavailability,
            time_vector,
            String[]  # No excluded feeders in this case
        )
    end
    
    return results
end

"""
    create_ac_network(config, time_vector::Vector) -> (_MSS.Network, Vector{Float64})

Create and solve AC network model.

# Arguments
- `config`: Configuration object
- `time_vector::Vector`: Time vector for analysis

# Returns
- `Tuple`: (AC Network object, AC network availability)
"""
function create_ac_network(config, time_vector::Vector)
    # Create AC grid STD
    ac_availability = config.simulation.ac_grid_availability * ones(length(time_vector))
    ac_std = _MSS.solvedSTD(
        prob = [ac_availability, 1 .- ac_availability],
        time = time_vector,
        power = [config.simulation.source_power, config.simulation.zero_power]
    )
    
    # Note: This is a simplified version. The full implementation would need
    # access to the AC/DC converter STD data from the original analysis
    @warn "AC network creation simplified - converter STD data needed for complete implementation"
    
    # Create network
    network = _MSS.Network()
    _MSS.add_user!(network, node = config.network.dc_bus_nodes)
    
    # Add AC sources
    for i in 1:config.network.ac_sources
        _MSS.add_sources!(network, node = i, std = ac_std)
    end
    
    # Solve network
    _MSS.solve!(network)
    
    # Extract availability (simplified - needs proper implementation with converter data)
    ac_network_availability = ac_availability  # Placeholder
    
    return network, ac_network_availability
end

"""
    print_analysis_summary(results::Dict)

Print summary of analysis results.

# Arguments
- `results::Dict`: Analysis results to summarize
"""
function print_analysis_summary(results::Dict)
    println("=== Bus Availability Analysis Summary ===")
    
    for (device_type, result) in results
        avg_availability = mean(result.availability_prob)
        min_availability = minimum(result.availability_prob)
        max_unavailability = maximum(result.unavailability_prob)
        
        println("$device_type:")
        println("  Average Availability: $(round(avg_availability, digits=6))")
        println("  Minimum Availability: $(round(min_availability, digits=6))")
        println("  Maximum Unavailability: $(round(max_unavailability, digits=8))")
        
        if !isempty(result.excluded_feeders)
            println("  Excluded Feeders: $(join(result.excluded_feeders, ", "))")
        end
        println()
    end
    
    println("==========================================")
end
