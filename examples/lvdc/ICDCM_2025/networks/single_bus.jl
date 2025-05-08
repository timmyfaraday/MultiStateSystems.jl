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

"""
    Load the required packages and modules for the simulation.
"""

include(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/functions/helper.jl"))

jld = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/single_bus_SemiMarkov.jld"), "r")

std_loads = jld["loads"]
std_sources = jld["sources"]

"""
    Create a source state transition diagram (STD) for the AC grid, create a dictionary of networks and solve each network to determine the instantaneous availabilities of the loads, sources and dc bus.
"""

AC_av = 0.99999

stdˢ = solvedSTD(prob = [ones(length(std_sources[1.0]["SSCB"]["S1"].props[:time]))*AC_av, ones(length(std_sources[1.0]["SSCB"]["S1"].props[:time]))*(1 .- AC_av)], time = collect(std_sources[1.0]["SSCB"]["S1"].props[:time]), power = [(Inf)u"MW", 0.0u"MW"]);

netw_1 = create_network(std_loads, std_sources, stdˢ)

solve_network_dict!(netw_1)

"""
    Find the network reliability of the dc bus and loads through the use of the `h` function, which is the transition frequency density function.
"""
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

h_load_2 = Dict(
    key_zone => Dict(
        key_type => Dict(
            key => sum(feeder.sprops[3][:h] .+ feeder.sprops[5][:h] for feeder in values(std_feeder)) .+
                   sum(source.sprops[3][:h] .+ source.sprops[5][:h] for source in values(std_sources[key_zone]["SSCB"])) .+
                   std_load.sprops[2][:h]
            for (key, std_load) in std_feeder
        )
        for (key_type, std_feeder) in std_CB
    )
    for (key_zone, std_CB) in std_loads
)


# # Create a DataFrame to store the h_load data
# df_h_load = DataFrame()

# # Populate the DataFrame with the data from h_load
# for (key_zone, h_zone) in h_load
#     for (key_type, h_cable) in h_zone
#         for (key, h_value) in h_cable
#             column_name = string(key_zone, "_", key_type, "_", key)
#             df_h_load[!, column_name] = h_value
#         end
#     end
# end

# Save the DataFrame to a CSV file
CSV.write(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/h_load_results.csv"), df_h_load)

function h_AB(h, P_zone, key_CB, key_feeder1, key_feeder2)
    fail1 = h[P_zone][key_CB][key_feeder1]
    fail2 = h_load[P_zone][key_CB][key_feeder2]

    U1 = 

    function extract_integer_from_string(s::String)
        match = match(r"\d+", s)
        return match !== nothing ? parse(Int, match.match) : nothing
    end

fail1 = h_load[0.999]["MCCB"]["C2"]
fail2 = h_load[0.999]["MCCB"]["C4"]

U1 = netw_1[0.999]["MCCB"].usr[3][:ugf].prb[1]
U2 = netw_1[0.999]["MCCB"].usr[5][:ugf].prb[1]

h_combined = DSP.conv(fail1*dt, U2) .+ DSP.conv(fail2*dt, U1)
dt = std_sources[1.0]["SSCB"]["S1"].props[:time][2]

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

h_combined = int1 .+ int2

CSV.write(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/h_combined_results.csv"), DataFrame(h_combined = h_combined[1:60:end]))

