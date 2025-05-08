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

include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/lvdc_netw_funcs.jl"))
include(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/functions/helper.jl"))

# Include the single bus SemiMarkov.jld in case the double bus network is not linked through cable in between the two zones.
jld_1 = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/1bus_SemiMarkov.jld"), "r")
# Include the double bus SemiMarkov.jld in case the double bus network is linked through cable in between the two zones.
jld_2 = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/2bus_SemiMarkov.jld"), "r")

std_loads_1 = jld_1["loads"]
std_sources_1 = jld_1["sources"]

std_loads_2 = jld_2["loads"]
std_sources_2 = jld_2["sources"]
std_bridge_2 = jld_2["bridge"] 

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

using DataFrames

# Create a DataFrame to store the data
function write_results(dict, red, input)
    for (key_zone, netw_CB) in dict
        for (key_CB, netw) in netw_CB
            for i in 1:length(netw.usr)            
                df[!, string(key_zone, "_", key_CB, "_User", i)] = netw.usr[i][:ugf].prb[input][1:red:end]
            end
        end
    end
end

df = DataFrame()
df[!, "Time"] = collect(std_sources_1[1.0]["SSCB"]["S1"].props[:time])[1:60:end]
write_results(double_bus_network_4, 60,1)

CSV.write(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/double_bus_network_results_far_double_SSCB_inverse.csv"), df)


h_bus_zone_1 = Dict()
for (key_zone, std_CB) in std_loads_1
    h_sum = Dict()
    for(key_type, std_feeder) in std_CB
        h_sum[key_type] = sum(feeder.sprops[3][:h] .+ feeder.sprops[5][:h] for (key,feeder) in std_feeder if key in Zone_1) .+
        .+ std_sources[key_zone]["SSCB"]["S1"].sprops[3][:h] .+ std_sources[key_zone]["SSCB"]["S1"].sprops[5][:h] 

    end
    h_bus_zone_1[key_zone] = h_sum
end

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


h_bus_zone_2 = Dict()
for (key_zone, std_CB) in std_loads_1
    h_sum = Dict()
    for(key_type, std_feeder) in std_CB
        h_sum[key_type] = sum(feeder.sprops[3][:h] .+ feeder.sprops[5][:h] for (key,feeder) in std_feeder if key in Zone_2) .+
        .+ std_sources[key_zone]["SSCB"]["S1"].sprops[3][:h] .+ std_sources[key_zone]["SSCB"]["S1"].sprops[5][:h] 

    end
    h_bus_zone_2[key_zone] = h_sum
end

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

fail1 = h_load_zone_1[0.999]["MCCB"]["C2"]
fail2 = h_load_zone_2[0.999]["MCCB"]["C4"]

U1 = double_bus_network_1[0.999]["MCCB"].usr[4][:ugf].prb[1]
U2 = double_bus_network_1[0.999]["MCCB"].usr[6][:ugf].prb[1]

int1 = zeros(length(fail1))
for i in 1:length(fail1)    
    for j in 1:i
        int1[i] += fail1[j] * dt * U2[length(fail1)-j+1] * MultiStateSystems.weights(i,j)
    end
end

int2 = zeros(length(fail1))
for i in 1:length(fail1)    
    for j in 1:i
        int1[i] += fail2[j] * dt * U1[length(fail1)-j+1] * MultiStateSystems.weights(i,j)
    end
end

h_combined_double_bus = int1 .+ int2


# Create a DataFrame to store the h_load data
df_h_load_1 = DataFrame()

# Populate the DataFrame with the data from h_load
for (key_zone, h_zone) in h_load_zone_1
    for (key_type, h_cable) in h_zone
        for (key, h_value) in h_cable
            column_name = string(key_zone, "_", key_type, "_", key)
            df_h_load_1[!, column_name] = h_value
        end
    end
end

# Save the DataFrame to a CSV file
CSV.write(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/h_load_zone1_results.csv"), df_h_load_1)

df_h_load_2 = DataFrame()
# Populate the DataFrame with the data from h_load
for (key_zone, h_zone) in h_load_zone_2
    for (key_type, h_cable) in h_zone
        for (key, h_value) in h_cable
            column_name = string(key_zone, "_", key_type, "_", key)
            df_h_load_2[!, column_name] = h_value
        end
    end
end
# Save the DataFrame to a CSV file
CSV.write(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/h_load_zone2_results.csv"), df_h_load_2)

# Save the h_combined_double_bus to a CSV file
CSV.write(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/h_combined_double_bus.csv"), DataFrame(h_combined_double_bus = h_combined_double_bus[1:60:end]))