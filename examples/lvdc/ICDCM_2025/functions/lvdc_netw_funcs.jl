################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/DC_Faults.jl"))

# Function to determine probability of repair before battery reserve time is reached
function battery_system_availability(i, T, ntw_av, μ, σ)
    # Calculate the probability that the battery does not run out of charge before
    # the power sources are repaired at time t. With a lognormal distribution for 
    # the repair time with mean μ and standard deviation σ.
    p_failure = 1 - ntw_av[i]
    p_repair_time_ge_T = _MSS.ccdf(LogNormal(μ, σ),T)
    return 1-p_failure*p_repair_time_ge_T
end

function calculate_bus_availability(exclude_keys::Union{Vector{String}, Nothing} = nothing)
    stdᵇᵘˢ = Dict()

    for (key, L_c) in L_tot
        stdᵇᵘˢ[key] = solvedSTD(
            prob = [1 .- sum([_MSS.get_sprop(value, :prob)[3] for (i, value) in std_s[key] if exclude_keys === nothing || !(key in exclude_keys)]),
                    sum([_MSS.get_sprop(value, :prob)[3] for (i, value) in std_s[key] if exclude_keys === nothing || !(key in exclude_keys)])],
            time = collect(time),
            power = [(Inf)u"MW", 0.0u"MW"]
        )
    end

    return stdᵇᵘˢ
end

function calculate_bus_availability_with_keys(include_keys::Vector{String})
    stdᵇᵘˢ = Dict()
    for (key, L_c) in L_tot
            stdᵇᵘˢ[key] = solvedSTD(
                prob = [1 .- sum([_MSS.get_sprop(std_s[key][i], :prob)[3] for i in include_keys]),
                        sum([_MSS.get_sprop(std_s[key][i], :prob)[3] for i in include_keys])],
                time = collect(time),
                power = [(Inf)u"MW", 0.0u"MW"]
            )
    end

    return stdᵇᵘˢ
end

function calculate_P(L_tot, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
# Calculate the probability of clearing a fault for each fault location and protection device
    P = Dict()
    for (device_type, feeder_lengths) in L_tot
        probabilities = Dict()

        if occursin("_", device_type)
            device_a, device_b = split(device_type, "_")
            for (feeder, lengths) in feeder_lengths
                prob_a = λ[device_a][2] == "CB" ? fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_a], λ[device_a][1]; return_mean = false) :
                                                fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_a][1]; return_mean = false)
                prob_b = λ[device_b][2] == "CB" ? fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_b], λ[device_b][1]; return_mean = false) :
                                                fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_b][1]; return_mean = false)
                probabilities[feeder] = mean(max.(prob_a, prob_b))
            end
        else
            for (feeder, lengths) in feeder_lengths
                probabilities[feeder] = λ[device_type][2] == "CB" ? fault_clear_prob(lengths, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_type], λ[device_type][1]) :
                                                                    fault_clear_prob_fuse(lengths, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_type][1])
            end
        end

        P[device_type] = probabilities
    end
    return P
end

function calculate_Pc(L_tot, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
    Pc = Dict()
    for (device_type, feeder_lengths) in L_tot
        probabilities = Dict()

        if occursin("_", device_type)
            device_a, device_b = split(device_type, "_")
            for (feeder, lengths) in feeder_lengths
                prob_a = λ[device_a][2] == "CB" ? fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_a], λ[device_a][1]; return_mean = false) :
                                                fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_a][1]; return_mean = false)
                prob_b = λ[device_b][2] == "CB" ? fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_b], λ[device_b][1]; return_mean = false) :
                                                fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_b][1]; return_mean = false)
                probabilities[feeder] = mean(max.(prob_a, prob_b))
            end
        else
            for (feeder, lengths) in feeder_lengths
                probabilities[feeder] = λ[device_type][2] == "CB" ? fault_clear_prob(maximum(lengths), L_p, C_b, V_DC, I_max, V_min, n, t_max, μ[device_type], λ[device_type][1]) :
                                                                    fault_clear_prob_fuse(maximum(lengths), L_p, C_b, V_DC, I²t, V_min, n, t_max, λ[device_type][1])
            end
        end

        Pc[device_type] = probabilities
    end
    return Pc
end