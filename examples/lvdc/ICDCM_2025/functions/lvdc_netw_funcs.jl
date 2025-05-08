################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using MultiStateSystems
using Unitful

const _MSS = MultiStateSystems

# Include additional functions for DC faults
include(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/functions/DC_Faults.jl"))

# """
#     calculate_bus_availability(exclude_keys::Union{Vector{String}, Nothing} = nothing)

# Calculate the availability of buses in the system, optionally excluding certain keys.

# # Arguments
# - `exclude_keys::Union{Vector{String}, Nothing}`: A vector of keys to exclude from the calculation, or `nothing` to include all keys.

# # Returns
# - A dictionary mapping feeder availabilities to their bus state transition diagrams (STDs).
# """
# function calculate_bus_availability(exclude_keys::Union{Vector{String}, Nothing} = nothing)
#     stdᵇᵘˢ = Dict()

#     for (key, L_c) in L_tot
#         stdᵇᵘˢ[key] = solvedSTD(
#             prob = [
#                 1 .- sum([_MSS.get_sprop(value, :prob)[3] for (i, value) in std_s[key] if exclude_keys === nothing || !(key in exclude_keys)]),
#                 sum([_MSS.get_sprop(value, :prob)[3] for (i, value) in std_s[key] if exclude_keys === nothing || !(key in exclude_keys)])
#             ],
#             time = collect(time),
#             power = [(Inf)u"MW", 0.0u"MW"]
#         )
#     end
#     return stdᵇᵘˢ
# end

# """
#     calculate_bus_availability_with_keys(include_keys::Vector{String})

# Calculate the availability of buses in the system, considering only the specified keys.

# # Arguments
# - `include_keys::Vector{String}`: A vector of keys to include in the calculation.

# # Returns
# - A dictionary mapping bus keys to their solved state transition diagrams (STDs).
# """
# function calculate_bus_availability_with_keys(include_keys::Vector{String})
#     stdᵇᵘˢ = Dict()

#     for (key, L_c) in L_tot
#         stdᵇᵘˢ[key] = solvedSTD(
#             prob = [
#                 1 .- sum([_MSS.get_sprop(std_s[key][i], :prob)[3] for i in include_keys]),
#                 sum([_MSS.get_sprop(std_s[key][i], :prob)[3] for i in include_keys])
#             ],
#             time = collect(time),
#             power = [(Inf)u"MW", 0.0u"MW"]
#         )
#     end
#     return stdᵇᵘˢ
# end

"""
    calculate_P(L_tot, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)

Calculate the probability of clearing a fault for each fault location and protection device.

# Arguments
- `L_tot`: Total lengths of feeders for each device type.
- `L_p`: Length of the protection zone.
- `C_b`: Capacitance of the bus.
- `V_DC`: DC voltage.
- `I_max`: Maximum current.
- `V_min`: Minimum voltage.
- `n`: Number of samples.
- `t_max`: Maximum time.
- `μ`: Mean time to clear faults for each device type.
- `λ`: Fault rate and type information for each device type.

# Returns
- A dictionary mapping device types to their fault-clearing probabilities for each feeder.
"""
function calculate_P(L_tot, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
    P = Dict()

    for (device_type, feeder_lengths) in L_tot
        probabilities = Dict()

        if occursin("_", device_type)
            device_a, device_b = split(device_type, "_")
            for (feeder, lengths) in feeder_lengths
                prob_a = λ[device_a][2] == "CB" ?
                    fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_a], λ[device_a][1]; return_mean = false) :
                    fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_a][1]; return_mean = false)

                prob_b = λ[device_b][2] == "CB" ?
                    fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_b], λ[device_b][1]; return_mean = false) :
                    fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_b][1]; return_mean = false)

                probabilities[feeder] = mean(max.(prob_a, prob_b))
            end
        else
            for (feeder, lengths) in feeder_lengths
                probabilities[feeder] = λ[device_type][2] == "CB" ?
                    fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_type], λ[device_type][1]) :
                    fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_type][1])
            end
        end
        P[device_type] = probabilities
    end
    return P
end

"""
    calculate_Pc(L_tot, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)

Calculate the probability of clearing the worst-case fault for each fault location and protection device.

# Arguments
- `L_tot`: Total lengths of feeders for each device type.
- `L_p`: Length of the protection zone.
- `C_b`: Capacitance of the bus.
- `V_DC`: DC voltage.
- `I_max`: Maximum current.
- `V_min`: Minimum voltage.
- `n`: Number of samples.
- `t_max`: Maximum time.
- `μ`: Mean time to clear faults for each device type.
- `λ`: Fault rate and type information for each device type.

# Returns
- A dictionary mapping device types to their worst-case fault-clearing probabilities for each feeder.
"""
function calculate_Pc(L_tot, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
    Pc = Dict()

    for (device_type, feeder_lengths) in L_tot
        probabilities = Dict()

        if occursin("_", device_type)
            device_a, device_b = split(device_type, "_")
            for (feeder, lengths) in feeder_lengths
                prob_a = λ[device_a][2] == "CB" ?
                    fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_a], λ[device_a][1]; return_mean = false) :
                    fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_a][1]; return_mean = false)

                prob_b = λ[device_b][2] == "CB" ?
                    fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_b], λ[device_b][1]; return_mean = false) :
                    fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_b][1]; return_mean = false)

                probabilities[feeder] = mean(max.(prob_a, prob_b))
            end
        else
            for (feeder, lengths) in feeder_lengths
                probabilities[feeder] = λ[device_type][2] == "CB" ?
                    fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_type], λ[device_type][1]) :
                    fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_type][1])
            end
        end
        Pc[device_type] = probabilities
    end
    return Pc
end

function create_network(std_loads, std_sources, std_ac_grid)
    single_bus_networks = Dict()

    # Iterate over the keys and values of the original dictionary   
    for (key_zone, std_CB) in std_loads
        single_bus_network_CB = Dict()
        for(key_type, std_feeder) in std_CB
            # Check if the value is a dictionary
            netw = Network()
            add_sources!(netw, node = [1,2], std = [std_ac_grid, std_ac_grid])
            add_components!(netw, edge = [(1,3), (2,3)],
                            std = [std_sources[key_zone]["SSCB"]["S1"], std_sources[key_zone]["SSCB"]["S2"]])
            # Calculate DC bus availability depending on the availability of the different feeders.
            Uᵦ = zeros(length(std_sources[1.0]["SSCB"]["S1"].props[:time]))
            for feeder in values(std_feeder)
                for state in values(feeder.sprops)
                    if occursin("U", state[:name])
                        Uᵦ += state[:prob]
                    end
                end
            end
            for feeder in values(std_sources[key_zone]["SSCB"])
                for state in values(feeder.sprops)
                    if occursin("U", state[:name])
                        Uᵦ += state[:prob]
                    end
                end
            end
            stdᵇ = solvedSTD(prob = [1 .- Uᵦ, Uᵦ], power = [(Inf)u"MW", 0.0u"MW"], time = collect(std_sources[1.0]["SSCB"]["S1"].props[:time]))
            add_component!(netw, node = 3, std = stdᵇ)
            add_components!(netw, edge = [(3, 3 + i) for i in 1:length(std_feeder)],
                            std = [std_l for std_l in values(std_feeder)])
            
            add_users!(netw, node = [3,4,5,6,7,8])  
            single_bus_network_CB[key_type] = netw      
        end
        single_bus_networks[key_zone] = single_bus_network_CB
    end
    return single_bus_networks
end

function create_split_network(std_loads, std_sources, std_ac_grid, Zone_1, Zone_2; std_bridge = nothing)
    double_bus_networks = Dict()
    # Iterate over the keys and values of the original dictionary   
    for (key_zone, std_CB) in std_loads
        double_bus_network_CB = Dict()
        for(key_type, std_feeder) in std_CB
            # Check if the value is a dictionary
            netw = Network()
            add_sources!(netw, node = [1,2], std = [std_ac_grid, std_ac_grid])
            add_components!(netw, edge = [(1,3), (2,4)],
                            std = [std_sources[key_zone]["SSCB"]["S1"], std_sources[key_zone]["SSCB"]["S2"]])
            # Calculate DC bus availability depending on the availability of the different feeders.
            Uᵦ₁ = sum(state[:prob] for (key, feeder) in std_feeder if key in Zone_1 for state in values(feeder.sprops) if occursin("U", state[:name])) .+
            sum(state[:prob] for state in values(std_sources[key_zone]["SSCB"]["S1"].sprops) if occursin("U", state[:name])) .+
            (std_bridge == nothing ? 
                key_zone * sum(state[:prob] for (key, feeder) in std_feeder if key in Zone_2 for state in values(feeder.sprops) if occursin("U3", state[:name])) :
                ("C" in Zone_1 ? sum(state[:prob] for state in values(std_bridge[key_zone]["SSCB"]["C"].sprops) if occursin("U", state[:name])) : (1 .- std_bridge[key_zone]["SSCB"]["C"].sprops[1][:prob])))
                
    
            Uᵦ₂ = sum(state[:prob] for (key, feeder) in std_feeder if key in Zone_2 for state in values(feeder.sprops) if occursin("U", state[:name])) .+
            sum(state[:prob] for state in values(std_sources[key_zone]["SSCB"]["S2"].sprops) if occursin("U", state[:name])) .+
            (std_bridge == nothing ? 
                key_zone * sum(state[:prob] for (key, feeder) in std_feeder if key in Zone_1 for state in values(feeder.sprops) if occursin("U3", state[:name])) :
                ("C" in Zone_2 ? sum(state[:prob] for state in values(std_bridge[key_zone]["SSCB"]["C"].sprops) if occursin("U", state[:name])) : (1 .- std_bridge[key_zone]["SSCB"]["C"].sprops[1][:prob])))    
            
            stdᵇ¹ = solvedSTD(prob = [1 .- Uᵦ₁, Uᵦ₁], power = [(Inf)u"MW", 0.0u"MW"], time = collect(std_sources[1.0]["SSCB"]["S1"].props[:time]))
            stdᵇ² = solvedSTD(prob = [1 .- Uᵦ₂, Uᵦ₂], power = [(Inf)u"MW", 0.0u"MW"], time = collect(std_sources[1.0]["SSCB"]["S1"].props[:time]))	
            add_component!(netw, node = 3, std = stdᵇ¹)
            add_component!(netw, node = 4, std = stdᵇ²)

            # create std of the SSCB connecting zone 1 and zone 2
            if std_bridge == nothing
                std_bridge_close = solvedSTD(prob = [ones(length(std_sources[1.0]["SSCB"]["S1"].props[:time]))], power = [(Inf)u"MW"], time = collect(std_sources[1.0]["SSCB"]["S1"].props[:time]))
                add_components!(netw, edge = [(3,4), (4,3)], 
                                  std = [std_bridge_close, std_bridge_close])
            else
                add_components!(netw, edge = [(3,4), (4,3)], 
                                  std = [std_bridge[key_zone]["SSCB"]["C"], std_bridge[key_zone]["SSCB"]["C"]])
            end
          

            add_components!(netw, edge = [(3, 5), (3,6)],
                            std = [std_feeder["C1"], std_feeder["C2"]])
            add_components!(netw, edge = [(4, 7), (4,8), (4,9)],
                            std = [std_feeder["C3"], std_feeder["C4"], std_feeder["C5"]])
                        
            add_users!(netw, node = [3,4,5,6,7,8,9])  
            double_bus_network_CB[key_type] = netw      
        end
        double_bus_networks[key_zone] = double_bus_network_CB
    end
    return double_bus_networks
end