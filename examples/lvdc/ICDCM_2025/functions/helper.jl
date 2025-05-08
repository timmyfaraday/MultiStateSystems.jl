################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using MultiStateSystems
using Unitful

"""
    create_reduced_std!(std::Dict)

Creates a reduced version of the given state transition dictionary (`std`).
The function iterates over the keys and values of the input dictionary, processes
each sub-dictionary, and generates a reduced state transition dictionary.

# Arguments
- `std::Dict`: The original state transition dictionary.

# Returns
- `Dict`: A reduced state transition dictionary.
"""
function create_reduced_std!(std::Dict)
    std_reduced = Dict()

    for (key, value) in std
        if isa(value, AbstractDict)
            reduced_std_i = Dict()
            for (cb, std_sol) in value
                reduced_std_i[cb] = fill_reduced_std(std_sol)
            end
            std_reduced[key] = reduced_std_i
        end
    end

    return std_reduced
end

"""
    fill_reduced_std(std_sol::Any)

Processes a state transition solution (`std_sol`) and generates a reduced state
transition dictionary with probabilities and hazard rates for each state.

# Arguments
- `std_sol::Any`: The state transition solution to process.

# Returns
- `solvedSTD`: A reduced state transition dictionary containing probabilities,
  hazard rates, and other properties for each state.
"""
function fill_reduced_std(std_sol::Any)
    time = std_sol.props[:time]
    num_states = 6
    state_probs = [zeros(length(time)) for _ in 1:num_states]
    state_hazards = [zeros(length(time)) * unit(std_sol.sprops[1][:h][1]) for _ in 1:num_states]

    for n in 1:length(std_sol.sprops)
        state_name = std_sol.sprops[n][:name]
        prob = std_sol.sprops[n][:prob]
        hazard = std_sol.sprops[n][:h]

        if occursin("A", state_name)
            state_probs[1] .+= prob
            state_hazards[1] .+= hazard
        elseif occursin("V1", state_name)
            state_probs[2] .+= vcat(0, prob[2:end])
            state_hazards[2] .+= hazard
        elseif occursin("U2", state_name)
            corrective_prob = state_conv(LogNormal(log(2.0)u"hr", 0.25u"hr"), hazard, time, 10000)
            state_probs[3] .+= corrective_prob
            state_hazards[3] .+= hazard
            state_probs[4] .+= vcat(0, prob[2:end]) .- corrective_prob
            state_hazards[4] .+= hazard
        elseif occursin("U3", state_name)
            corrective_prob = state_conv(LogNormal(log(20.0)u"hr", 0.05u"hr"), hazard, time, 10000)
            state_probs[5] .+= corrective_prob
            state_hazards[5] .+= hazard
            state_probs[6] .+= vcat(0, prob[2:end]) .- corrective_prob
            state_hazards[6] .+= vcat(0 * unit(std_sol.sprops[1][:h][1]), hazard[1:end-1])
        end
    end

    return solvedSTD(
        prob = state_probs,
        time = collect(time),
        power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW"],
        name = ["A", "V1", "U2", "V2", "U3", "V3"],
        h = state_hazards
    )
end

# function process_nested_dict(input_dict, func)
#     # Recursive function to process the nested dictionary
#     function process_recursive(dict)
#         result = Dict()
#         for (key, value) in dict
#             if isa(value, AbstractDict)
#                 # If the value is a dictionary, recurse
#                 result[key] = process_recursive(value)
#             else
#                 # Apply the function to the lowest level value
#                 result[key] = func(value)
#             end
#         end
#         return result
#     end

#     # Call the recursive function on the input dictionary
#     return process_recursive(input_dict)
# end

# using SMTPClient

# function send_email_notification()
#     opt = SMTPClient.SendOptions(
#             isSSL = true,
#             username = "emmersglenn@gmail.com",
#             passwd = "kysn tqhv ubyv cikr")   
        
#     body = IOBuffer(
#     "From: You <emmersglenn@gmail.com>\r\n" *
#     "To: glenn.emmers@kuleuven.be\r\n" *
#     "Subject: Julia Finish\r\n" *
#     "\r\n" *
#     "Your code is done\r\n")
#     url = "smtps://smtp.gmail.com:465"
#     rcpt = ["<glenn.emmers@kuleuven.be>"]
#     from = "<emmersglenn@gmail.com>"
#     return SMTPClient.send(url, rcpt, from, body, opt)
# end

"""
    solve_network_dict!(netw_1::AbstractDict)

Recursively solves a nested dictionary of network components. This function traverses 
the input dictionary `netw_1`, and for each key-value pair, it checks if the value is 
another dictionary. If so, it recursively calls itself on that dictionary. Otherwise, 
it assumes the value is a solvable object and calls `solve!` on it.

# Arguments
- `netw_1::AbstractDict`: A dictionary where values can either be nested dictionaries 
  or objects that implement the `solve!` method.

# Behavior
- If a value in the dictionary is another dictionary, the function is called recursively 
  on that value.
- If a value is not a dictionary, the function assumes it is a solvable object and 
  applies the `solve!` function to it.

"""
function solve_network_dict!(netw_1)
    for (key, value) in netw_1
        if isa(value, AbstractDict)
            solve_network_dict!(value)
        else
            solve!(value)
        end
    end
end