using MultiStateSystems
using Unitful

function create_reduced_std!(std)
    # Create a new dictionary to store the reduced std
    std_reduced = Dict()

    # Iterate over the keys and values of the original dictionary
    for (key, value) in std
        # Check if the value is a dictionary
        if isa(value, AbstractDict)
            # Create a new dictionary for the reduced std of the current key
            reduced_std_i = Dict()
            for (cb, std_sol) in value
                # Create a new STD with the same properties as the original one
                reduced_std_i[cb] = fill_reduced_std(std_sol)
            end
            # Store the reduced std in the new dictionary
            std_reduced[key] = reduced_std_i
        end
    end

    return std_reduced
end
    
function fill_reduced_std(std_sol)
    time = std_sol.props[:time]
    state_1_prob = zeros(length(time))
    state_2_prob = zeros(length(time))
    state_3_prob = zeros(length(time)) 
    state_4_prob = zeros(length(time))
    state_5_prob = zeros(length(time))
    state_6_prob = zeros(length(time))
    for n in 1:length(std_sol.sprops)
        if occursin("A", std_sol.sprops[n][:name])
            state_1_prob .+= std_sol.sprops[n][:prob]
        elseif occursin("V1", std_sol.sprops[n][:name])
            state_2_prob .+= vcat(0, std_sol.sprops[n][:prob][2:end])
        elseif occursin("U2", std_sol.sprops[n][:name])
            corrective_state_prob_1 = state_conv(LogNormal(log(2.0)u"hr", 0.25u"hr"), std_sol.sprops[n][:h], time, 10000)
            state_3_prob .+= corrective_state_prob_1
            state_4_prob .+= vcat(0, std_sol.sprops[n][:prob][2:end]).-corrective_state_prob_1

        elseif occursin("U3", std_sol.sprops[n][:name])
            corrective_state_prob_2 = state_conv(LogNormal(log(2.0)u"hr", 0.25u"hr"), std_sol.sprops[n][:h], time, 10000)
            state_5_prob .+= corrective_state_prob_2
            state_6_prob .+= vcat(0, std_sol.sprops[n][:prob][2:end]).-corrective_state_prob_2
        end
    end

    return solvedSTD(prob = [state_1_prob, state_2_prob, state_3_prob, state_4_prob, state_5_prob, state_6_prob], time = collect(time), power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW"], name  = ["A", "U1", "U2","U3", "V2","V3"])
end

function process_nested_dict(input_dict, func)
    # Recursive function to process the nested dictionary
    function process_recursive(dict)
        result = Dict()
        for (key, value) in dict
            if isa(value, AbstractDict)
                # If the value is a dictionary, recurse
                result[key] = process_recursive(value)
            else
                # Apply the function to the lowest level value
                result[key] = func(value)
            end
        end
        return result
    end

    # Call the recursive function on the input dictionary
    return process_recursive(input_dict)
end

using SMTPClient

function send_email_notification()
    opt = SMTPClient.SendOptions(
            isSSL = true,
            username = "emmersglenn@gmail.com",
            passwd = "kysn tqhv ubyv cikr")   
        
    body = IOBuffer(
    "From: You <emmersglenn@gmail.com>\r\n" *
    "To: glenn.emmers@kuleuven.be\r\n" *
    "Subject: Julia Finish\r\n" *
    "\r\n" *
    "Your code is done\r\n")
    url = "smtps://smtp.gmail.com:465"
    rcpt = ["<glenn.emmers@kuleuven.be>"]
    from = "<emmersglenn@gmail.com>"
    return SMTPClient.send(url, rcpt, from, body, opt)
end