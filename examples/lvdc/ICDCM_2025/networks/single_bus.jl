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
include(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/functions/helper.jl"))

jld = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/single_bus_SemiMarkov.jld"), "r")

std_loads = jld["loads"]
std_sources = jld["sources"]

## Create a source state transition diagram (STD) for the AC grid, create a dictionary of networks and solve each network to determine the instantaneous availabilities of the loads, sources and dc bus.
AC_av = 0.99999

# State transition diagram of the AC grid as a source
stdˢ = solvedSTD(prob = [ones(length(std_sources[1.0]["SSCB"]["S1"].props[:time]))*AC_av, ones(length(std_sources[1.0]["SSCB"]["S1"].props[:time]))*(1 .- AC_av)], time = collect(std_sources[1.0]["SSCB"]["S1"].props[:time]), power = [(Inf)u"MW", 0.0u"MW"]);

netw_1 = create_network(std_loads, std_sources, stdˢ)

solve_network_dict!(netw_1)

## Find the network reliability of the dc bus and loads through the use of the `h` function, which is the transition frequency density function.

## Determine the reliability of every bus in the dictionary of networks
h_bus = Dict()
for (key_zone, std_CB) in std_loads
    h_sum = Dict()
    for(key_type, std_feeder) in std_CB
        h_sum[key_type] = sum(feeder.sprops[3][:h] .+ feeder.sprops[5][:h] for feeder in values(std_feeder))
        .+ sum(source.sprops[3][:h] .+ source.sprops[5][:h] for source in values(std_sources[key_zone]["SSCB"]))
    end
    h_bus[key_zone] = h_sum
end

## Determine the reliability of every load connected to every feeder in the dictionary of networks
h_load = Dict()
for (key_zone, std_CB) in std_loads
    h_zone = Dict()
    for (key_type, std_feeder) in std_CB
        h_cable = Dict()
        for (key, std_load) in std_feeder
            h_cable[key] = h_bus[key_zone][key_type] .+ std_load.sprops[2][:h] 
        end
        h_zone[key_type] = h_cable
    end
    h_load[key_zone] = h_zone
end

"""
    h_AB(h, P_zone, key_CB, key_feeder1, key_feeder2, user1, user2)

Compute the combined effect of failures and user probabilities in a power network.

# Arguments
- `h::Dict`: A nested dictionary containing failure data for different zones, circuit breakers, and feeders.
- `P_zone::Symbol`: The key representing the power zone in the `h` dictionary.
- `key_CB::Symbol`: The key representing the circuit breaker in the `h` dictionary.
- `key_feeder1::Symbol`: The key representing the first feeder in the `h` dictionary.
- `key_feeder2::Symbol`: The key representing the second feeder in the `h` dictionary.
- `user1::Symbol`: The key representing the first user in the network.
- `user2::Symbol`: The key representing the second user in the network.

# Returns
- `Vector{Float64}`: A vector representing the combined integral of failure probabilities weighted by user probabilities and system weights.

# Details
The function calculates two integrals:
1. The effect of failures in `key_feeder1` weighted by the probability of `user2`.
2. The effect of failures in `key_feeder2` weighted by the probability of `user1`.

These integrals are computed using a nested loop over the failure data and are combined element-wise to produce the final result.

# Notes
- The function assumes the existence of global variables `h_load`, `netw_1`, `dt`, and `_MSS.weights`.
- The `_MSS.weights(i, j)` function is used to compute weights for the integrals.
"""
function h_AB_single_bus(h, P_zone, key_CB, key_feeder1, key_feeder2, user1, user2)
    fail1 = h[P_zone][key_CB][key_feeder1]
    fail2 = h_load[P_zone][key_CB][key_feeder2]

    U1 = netw_1[P_zone][key_CB].usr[user1][:ugf].prb[1]
    U2 = netw_1[P_zone][key_CB].usr[user2][:ugf].prb[1]

    int1 = zeros(length(fail1))
    for i in 1:length(fail1)    
        for j in 1:i
            int1[i] += fail1[j] * dt * U2[length(fail1)-j+1] * _MSS.weights(i,j)
        end
    end

    int2 = zeros(length(fail1))
    for i in 1:length(fail1)    
        for j in 1:i
            int2[i] += fail2[j] * dt * U1[length(fail1)-j+1] * _MSS.weights(i,j)
        end
    end

    return int1 .+ int2
end

h_combined = h_AB_single_bus(h_load, P_zone, "MCCB", "C2", "C4", 3, 5)


