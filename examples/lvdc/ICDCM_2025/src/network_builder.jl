################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Network Builder Module

This module provides functions for constructing LVDC network models using the 
MultiStateSystems.jl framework. It includes functions for building single bus
and double bus network configurations.

## Main Functions
- `build_single_bus_network`: Construct single bus network topology
- `build_double_bus_network`: Construct double bus network topology
- `create_ac_grid_std`: Create AC grid source STD
- `calculate_bus_availability`: Calculate DC bus availability from component STDs
"""

using MultiStateSystems
using Unitful

"""
    create_ac_grid_std(time_vector, availability::Float64=0.99999)

Create a state transition diagram for the AC grid source.

# Arguments
- `time_vector`: Time vector for the STD
- `availability::Float64`: AC grid availability (default: 0.99999)

# Returns
- `solvedSTD`: AC grid STD with specified availability
"""
function create_ac_grid_std(time_vector, availability::Float64=0.99999)
    return solvedSTD(
        prob = [ones(length(time_vector)) * availability, 
                ones(length(time_vector)) * (1 - availability)],
        time = collect(time_vector),
        power = [(Inf)u"MW", 0.0u"MW"]
    )
end

"""
    calculate_bus_availability(component_stds, zone_feeders=nothing)

Calculate DC bus availability based on component STDs.

# Arguments
- `component_stds`: Dictionary of component STDs
- `zone_feeders`: Optional list of feeders to include (default: all)

# Returns
- `solvedSTD`: DC bus STD with calculated availability
"""
function calculate_bus_availability(component_stds, zone_feeders=nothing)
    # Get time vector from first component
    first_std = first(values(component_stds))
    time_vector = first_std.props[:time]
    
    # Calculate total unavailability
    total_unavailability = zeros(length(time_vector))
    
    for (feeder_id, std) in component_stds
        if zone_feeders === nothing || feeder_id in zone_feeders
            # Sum unavailability states (states with "U" in name)
            for state in values(std.sprops)
                if occursin("U", state[:name])
                    total_unavailability .+= state[:prob]
                end
            end
        end
    end
    
    # Bus availability is complement of total unavailability
    availability = 1 .- total_unavailability
    
    return solvedSTD(
        prob = [availability, total_unavailability],
        time = collect(time_vector),
        power = [(Inf)u"MW", 0.0u"MW"]
    )
end

"""
    build_single_bus_network(load_stds, source_stds, ac_grid_std)

Build a single bus LVDC network topology.

# Arguments
- `load_stds`: Dictionary of load feeder STDs
- `source_stds`: Dictionary of source feeder STDs  
- `ac_grid_std`: AC grid STD

# Returns
- `Network`: Configured single bus network
"""
function build_single_bus_network(load_stds, source_stds, ac_grid_std)
    network = Network()
    
    # Add AC sources
    add_sources!(network, node=[1, 2], std=[ac_grid_std, ac_grid_std])
    
    # Add source converters (AC to DC)
    source_list = collect(values(source_stds))
    add_components!(network, 
                   edge=[(1, 3), (2, 3)], 
                   std=source_list)
    
    # Calculate DC bus availability
    all_component_stds = merge(load_stds, source_stds)
    bus_std = calculate_bus_availability(all_component_stds)
    add_component!(network, node=3, std=bus_std)
    
    # Add load feeders
    load_list = collect(values(load_stds))
    num_loads = length(load_list)
    add_components!(network,
                   edge=[(3, 3+i) for i in 1:num_loads],
                   std=load_list)
    
    # Add users (load connection points)
    user_nodes = [3; 3 .+ (1:num_loads)]
    add_users!(network, node=user_nodes)
    
    return network
end

"""
    build_double_bus_network(load_stds, source_stds, ac_grid_std, zone_config; bridge_stds=nothing)

Build a double bus LVDC network topology.

# Arguments
- `load_stds`: Dictionary of load feeder STDs
- `source_stds`: Dictionary of source feeder STDs
- `ac_grid_std`: AC grid STD
- `zone_config`: Dictionary defining which feeders belong to which zone
- `bridge_stds`: Optional bridge connection STDs

# Returns
- `Network`: Configured double bus network
"""
function build_double_bus_network(load_stds, source_stds, ac_grid_std, zone_config; bridge_stds=nothing)
    network = Network()
    zone_1_feeders = zone_config["Zone_1"]
    zone_2_feeders = zone_config["Zone_2"]
    
    # Add AC sources
    add_sources!(network, node=[1, 2], std=[ac_grid_std, ac_grid_std])
    
    # Add source converters to separate buses
    source_list = collect(values(source_stds))
    add_components!(network,
                   edge=[(1, 3), (2, 4)],
                   std=source_list)
    
    # Calculate bus availabilities for each zone
    zone_1_stds = Dict(k => v for (k, v) in load_stds if k in zone_1_feeders)
    zone_1_stds = merge(zone_1_stds, Dict("S1" => source_stds["S1"]))
    if bridge_stds !== nothing && "C" in zone_1_feeders
        zone_1_stds["C"] = bridge_stds["C"]
    end
    
    zone_2_stds = Dict(k => v for (k, v) in load_stds if k in zone_2_feeders)
    zone_2_stds = merge(zone_2_stds, Dict("S2" => source_stds["S2"]))
    if bridge_stds !== nothing && "C" in zone_2_feeders
        zone_2_stds["C"] = bridge_stds["C"]
    end
    
    # Handle cross-zone dependencies if no physical bridge
    if bridge_stds === nothing
        # Add cross-zone fault propagation
        for P_zone in [0.999, 1.0]  # This should come from zone_config
            zone_1_unavail = calculate_zone_unavailability(zone_1_stds, zone_2_stds, P_zone)
            zone_2_unavail = calculate_zone_unavailability(zone_2_stds, zone_1_stds, P_zone)
        end
    end
    
    bus_1_std = calculate_bus_availability(zone_1_stds)
    bus_2_std = calculate_bus_availability(zone_2_stds)
    
    add_component!(network, node=3, std=bus_1_std)
    add_component!(network, node=4, std=bus_2_std)
    
    # Add bridge connection
    if bridge_stds !== nothing
        bridge_std = bridge_stds["SSCB"]["C"]
        add_components!(network,
                       edge=[(3, 4), (4, 3)],
                       std=[bridge_std, bridge_std])
    else
        # Create ideal bridge (always closed)
        time_vector = bus_1_std.props[:time]
        ideal_bridge = solvedSTD(
            prob = [ones(length(time_vector))],
            time = collect(time_vector),
            power = [(Inf)u"MW"]
        )
        add_components!(network,
                       edge=[(3, 4), (4, 3)],
                       std=[ideal_bridge, ideal_bridge])
    end
    
    # Add load feeders to appropriate buses
    node_counter = 5
    for feeder_id in zone_1_feeders
        if haskey(load_stds, feeder_id)
            add_component!(network,
                          edge=(3, node_counter),
                          std=load_stds[feeder_id])
            node_counter += 1
        end
    end
    
    for feeder_id in zone_2_feeders
        if haskey(load_stds, feeder_id)
            add_component!(network,
                          edge=(4, node_counter),
                          std=load_stds[feeder_id])
            node_counter += 1
        end
    end
    
    # Add users
    user_nodes = vcat(3, 4, collect(5:node_counter-1)...)
    add_users!(network, node=user_nodes)
    
    return network
end

"""
    calculate_zone_unavailability(primary_zone_stds, secondary_zone_stds, P_zone)

Calculate zone unavailability considering cross-zone dependencies.

# Arguments
- `primary_zone_stds`: STDs for the primary zone
- `secondary_zone_stds`: STDs for the secondary zone
- `P_zone`: Zone protection factor

# Returns
- `Vector`: Zone unavailability over time
"""
function calculate_zone_unavailability(primary_zone_stds, secondary_zone_stds, P_zone)
    # Get time vector
    first_std = first(values(primary_zone_stds))
    time_vector = first_std.props[:time]
    
    # Calculate primary zone unavailability
    primary_unavail = zeros(length(time_vector))
    for std in values(primary_zone_stds)
        for state in values(std.sprops)
            if occursin("U", state[:name])
                primary_unavail .+= state[:prob]
            end
        end
    end
    
    # Add cross-zone fault propagation
    secondary_propagation = zeros(length(time_vector))
    for std in values(secondary_zone_stds)
        for state in values(std.sprops)
            if occursin("U3", state[:name])  # Only severe faults propagate
                secondary_propagation .+= state[:prob] * P_zone
            end
        end
    end
    
    return primary_unavail .+ secondary_propagation
end

"""
    solve_network!(network)

Solve a network to calculate user availability.

# Arguments
- `network`: Network object to solve

# Modifies
- `network`: Network is solved in-place
"""
function solve_network!(network)
    _MSS.solve!(network)
end

"""
    solve_multiple_networks!(networks)

Solve multiple networks in a nested dictionary structure.

# Arguments
- `networks`: Nested dictionary of networks

# Modifies
- `networks`: All networks are solved in-place
"""
function solve_multiple_networks!(networks)
    for (key, value) in networks
        if isa(value, AbstractDict)
            solve_multiple_networks!(value)
        else
            solve_network!(value)
        end
    end
end
