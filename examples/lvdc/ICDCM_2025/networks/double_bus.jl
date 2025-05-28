################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using MultiStateSystems
using Unitful
using JLD2

const _MSS = MultiStateSystems

## Load the required packages and modules for the simulation.
include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/lvdc_netw_funcs.jl"))
include(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/functions/helper.jl"))

# Include the single bus SemiMarkov.jld in case the double bus network is not linked through cable in between the two zones.
jld_1 = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/single_bus_SemiMarkov.jld"), "r")
# Include the double bus SemiMarkov.jld in case the double bus network is linked through cable in between the two zones.
jld_2 = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/double_bus_SemiMarkov.jld"), "r")

std_loads_1 = jld_1["loads"]
std_sources_1 = jld_1["sources"]

std_loads_2 = jld_2["loads"]
std_sources_2 = jld_2["sources"]
std_bridge_2 = jld_2["bridge"] 

## Create a source state transition diagram (STD) for the AC grid, create a dictionary of networks and solve each network to determine the instantaneous availabilities of the loads, sources and dc bus.
AC_av = 0.99999

# State transition diagram of the AC grid as a source
stdˢ = solvedSTD(prob = [ones(length(std_sources_1[1.0]["SSCB"]["S1"].props[:time]))*AC_av, ones(length(std_sources_1[1.0]["SSCB"]["S1"].props[:time]))*(1 .- AC_av)], time = collect(std_sources_1[1.0]["SSCB"]["S1"].props[:time]), power = [(Inf)u"MW", 0.0u"MW"]);

# Feeders connected to zone 1
Zone_1 = ["C1", "C2"]
# Feeders connected to zone 2
Zone_2 = ["C3", "C4", "C5"]

double_bus_network_1 = create_split_network(std_loads_1, std_sources_1, stdˢ, Zone_1, Zone_2; std_bridge = nothing)

double_bus_network_2 = create_split_network(std_loads_2, std_sources_2, stdˢ, ["C1", "C2", "C"], ["C3", "C4", "C5"]; std_bridge = std_bridge_2)

double_bus_network_3 = create_split_network(std_loads_2, std_sources_2, stdˢ, ["C1", "C2"], ["C3", "C4", "C5", "C"]; std_bridge = std_bridge_2)

double_bus_network_4 = create_split_network(std_loads_2, std_sources_2, stdˢ, ["C1", "C2", "C"], ["C3", "C4", "C5", "C"]; std_bridge = std_bridge_2)

solve_network_dict!(double_bus_network_1)
solve_network_dict!(double_bus_network_2)
solve_network_dict!(double_bus_network_3)
solve_network_dict!(double_bus_network_4)

## Find the network reliability of the dc bus and loads through the use of the `h` function, which is the transition frequency density function.

## Determine the reliability of every zone 1 in the dictionary of networks
h_bus_zone_1 = Dict()
for (key_zone, std_CB) in std_loads_1
    h_sum = Dict()
    for(key_type, std_feeder) in std_CB
        h_sum[key_type] = sum(feeder.sprops[3][:h] .+ feeder.sprops[5][:h] for (key,feeder) in std_feeder if key in Zone_1) .+
        .+ std_sources[key_zone]["SSCB"]["S1"].sprops[3][:h] .+ std_sources[key_zone]["SSCB"]["S1"].sprops[5][:h] 

    end
    h_bus_zone_1[key_zone] = h_sum
end

## Determine the reliability of every load connected to every feeder in zone 1 in the dictionary of networks
h_load_zone_1 = Dict()
for (key_zone, std_CB) in std_loads_1
    h_zone = Dict()
    for (key_type, std_feeder) in std_CB
        h_cable = Dict()
        for (key, std_load) in std_feeder; key in Zone_1
            h_cable[key] = h_bus_zone_1[key_zone][key_type] .+ std_load.sprops[2][:h] 
        end
        h_zone[key_type] = h_cable
    end
    h_load_zone_1[key_zone] = h_zone
end

## Determine the reliability of every zone 2 in the dictionary of networks
h_bus_zone_2 = Dict()
for (key_zone, std_CB) in std_loads_1
    h_sum = Dict()
    for(key_type, std_feeder) in std_CB
        h_sum[key_type] = sum(feeder.sprops[3][:h] .+ feeder.sprops[5][:h] for (key,feeder) in std_feeder if key in Zone_2) .+
        .+ std_sources[key_zone]["SSCB"]["S1"].sprops[3][:h] .+ std_sources[key_zone]["SSCB"]["S1"].sprops[5][:h] 

    end
    h_bus_zone_2[key_zone] = h_sum
end

## Determine the reliability of every load connected to every feeder in zone 2 in the dictionary of network
h_load_zone_2 = Dict()
for (key_zone, std_CB) in std_loads_1
    h_zone = Dict()
    for (key_type, std_feeder) in std_CB
        h_cable = Dict()
        for (key, std_load) in std_feeder; key in Zone_2
            h_cable[key] = h_bus_zone_2[key_zone][key_type] .+ std_load.sprops[2][:h] 
        end
        h_zone[key_type] = h_cable
    end
    h_load_zone_2[key_zone] = h_zone
end

"""
    h_AB_double_bus(h_zone1::Dict, h_zone2::Dict, network::Dict, key_zone, key_CB::Symbol, key_feeder1::Symbol, key_feeder2::Symbol, user1::Symbol, user2::Symbol) -> Vector

Compute the combined failure probabilities for a double-bus network configuration.

# Arguments
- `h_zone1::Dict`: A dictionary containing failure data for the first zone.
- `h_zone2::Dict`: A dictionary containing failure data for the second zone.
- `network::Dict`: A dictionary representing the network structure and user data.
- `key_zone`: The key identifying the zone in the dictionaries.
- `key_CB::Symbol`: The key identifying the circuit breaker in the zone.
- `key_feeder1::Symbol`: The key identifying the feeder in the first zone.
- `key_feeder2::Symbol`: The key identifying the feeder in the second zone.
- `user1::Symbol`: The key identifying the first user in the network.
- `user2::Symbol`: The key identifying the second user in the network.

# Returns
- `Vector`: A vector containing the combined failure probabilities for the double-bus network.

# Details
This function calculates the combined failure probabilities for a double-bus network by integrating failure data from two zones. It uses user probabilities (`U1` and `U2`) and failure data (`fail1` and `fail2`) to compute two integrals (`int1` and `int2`). These integrals are then combined element-wise to produce the final result.

The integration process involves iterating over the failure data and applying weights from `_MSS.weights` to compute the contributions from each failure event.

# Notes
- The variable `dt` is assumed to be defined globally and represents the time step for integration.
- `_MSS.weights` is a function or data structure that provides weights for the integration process.
"""
function h_AB_double_bus(h_zone1::Dict, h_zone2::Dict, network::Dict, key_zone, key_CB::Symbol, key_feeder1::Symbol, key_feeder2::Symbol, user1::Symbol, user2::Symbol)
    # Extract the failure data for the specified zone and circuit breaker
    fail1 = h_zone1[key_zone][key_CB][key_feeder1]
    fail2 = h_zone2[key_zone][key_CB][key_feeder2]

    # Extract the user probabilities for the specified users
    U1 = network[key_zone][key_CB].usr[user1][:ugf].prb[1]
    U2 = network[key_zone][key_CB].usr[user2][:ugf].prb[1]

    # Initialize the integrals
    int1 = zeros(length(fail1))
    int2 = zeros(length(fail2))

    # Calculate the first integral
    for i in 1:length(fail1)    
        for j in 1:i
            int1[i] += fail1[j] * dt * U2[length(fail1)-j+1] * _MSS.weights(i,j)
        end
    end

    # Calculate the second integral
    for i in 1:length(fail2)    
        for j in 1:i
            int2[i] += fail2[j] * dt * U1[length(fail2)-j+1] * _MSS.weights(i,j)
        end
    end

    # Combine the two integrals element-wise
    return int1 .+ int2
end

h_combined_double_bus = h_AB_double_bus(h_load_zone_1, h_load_zone_2, double_bus_network_1, 1.0, "SSCB", "C2", "C4", 4, 6)