################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

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
        if key in include_keys
            stdᵇᵘˢ[key] = solvedSTD(
                prob = [1 .- sum([_MSS.get_sprop(std_s[key][i], :prob)[3] for i in include_keys]),
                        sum([_MSS.get_sprop(std_s[key][i], :prob)[3] for i in include_keys])],
                time = collect(time),
                power = [(Inf)u"MW", 0.0u"MW"]
            )
        end
    end

    return stdᵇᵘˢ
end