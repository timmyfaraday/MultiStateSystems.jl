################################################################################
#  Copyright 2023, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
This document determines the fault current of an LVDC branch, depending on 
fault location and time.
The system is simplified towards an RLC circuit, which looks as follows:

|------R-------L-------|
|                      |
C                      |<-- Short circuit
|                      |
|------R-------L-------|
"""
# Load pkgs
using MultiStateSystems
using Unitful
using Statistics
using QuadGK
using Roots

"""
    sc_values(L_vector::Union{AbstractVector, Int}, L_p::Float64, C_b::Float64, V_DC::Int64) -> Vector{Function}

Compute the short-circuit current functions for a given set of cable lengths.

# Arguments
- `L_vector::Union{AbstractVector, Int}`: A vector of cable lengths (in meters) or a single integer representing a single cable length.
- `L_p::Float64`: The inductance of the power source (in Henries).
- `C_b::Float64`: The capacitance of the DC bus (in Farads).
- `V_DC::Int64`: The DC voltage (in Volts).

# Returns
- `Vector{Function}`: A vector of functions, where each function represents the short-circuit current as a function of time `t` (in seconds) for the corresponding cable length in `L_vector`.

# Details
The function calculates the short-circuit current for each cable length in `L_vector` based on the electrical parameters provided. It distinguishes between two cases:
1. **Underdamped response**: Occurs when the resistance squared (`R^2`) is less than `4 * L_total / C_b`. The current is modeled as a damped sinusoidal function.
2. **Overdamped response**: Occurs otherwise. The current is modeled as a combination of exponential terms.

To avoid numerical overflow in the overdamped case, the function includes a safeguard for large values of `ω * t`.

# Notes
- The cable parameters (inductance per unit length `l` and resistance per unit length `r`) are hardcoded as `2.94e-7 H/m` and `3.08e-3 Ω/m`, respectively.
- If `L_vector` is provided as an integer, it is converted into a single-element vector internally.
"""
function sc_values(L_vector::Union{AbstractVector, Int}, L_p::Float64, C_b::Float64, V_DC::Int64)
    # Cable parameters
    l = 2.94e-7  # H/m  
    r = 3.08e-3  # Ω/m

    # Handle the case where L_vector is an integer
    if isa(L_vector, Int)
        L_vector = [L_vector]
    end

    # Array to store the functions
    i_functions = Vector{Function}(undef, length(L_vector))

    for (idx, L) in enumerate(L_vector)
        R = 2 * r * L  # Resistance
        L_total = 2 * l * L + L_p
        if R^2 < 4 * L_total / C_b  # underdamped
            ω = sqrt(1 / (L_total * C_b) - (R / (2 * L_total))^2)
            i = t -> V_DC / (ω * L_total) * exp(-R / (2 * L_total) * t) * sin(ω * t)
        else  # overdamped
            ω = sqrt((R / (2 * L_total))^2 - 1 / (L_total * C_b))
            i = t -> begin
                exp_term = exp(-R / (2 * L_total) * t)
                sinh_term = exp(ω * t) - exp(-ω * t)
                if ω * t > 700
                    sinh_term = 2 * exp(ω * t - 700) * exp(700)  # Avoid overflow
                end
                V_DC / (2 * ω * L_total) * exp_term * sinh_term
            end
        end
        i_functions[idx] = i
    end
    return i_functions
end

"""
    Vmin(i_functions, C_b, t_int, V_DC, V_min)

Calculates the time at which the voltage across a DC bus drops below a specified minimum value (`V_min`) for a given set of current functions.

# Arguments
- `i_functions::Vector{Function}`: A vector of current functions, where each function represents the current as a function of time.
- `C_b::Float64`: The capacitance of the DC bus in Farads.
- `t_int::Unitful.Quantity`: The time interval over which the calculation is performed, with units of time.
- `V_DC::Float64`: The initial DC voltage in Volts.
- `V_min::Float64`: The minimum voltage threshold in Volts.

# Returns
- `Vector{Union{Float64, Nothing}}`: A vector containing the time at which the voltage drops below `V_min` for each current function. If no such time is found for a given function, the corresponding entry will be `Inf`.

# Notes
- The voltage function is computed as `v(t) = V_DC - ∫(i(t) / C_b) dt` over the interval `[0, t]`.
- The `find_zero` function is used with the Bisection method to find the root of the equation `v(t) - V_min = 0`.
- If no root is found for a given current function, the result for that function is set to `Inf`.
- If the root-finding process fails (e.g., due to numerical issues), the function attempts to estimate the time using a linear regression based on previous valid results. If no valid results exist, the time is set to `Inf`.
- The search interval is dynamically adjusted if the computed time approaches the upper bound of the interval.

# Exceptions
- If an error occurs during root finding for a specific current function, it is caught and handled by assigning `Inf` to the corresponding entry in the result vector.
"""

function Vmin(i_functions, C_b, t_int, V_DC, V_min)
    t_vmin = Vector{Union{Float64, Nothing}}(undef, length(i_functions))
    interval = (0.0, ustrip(t_int))

    for (idx, i) in enumerate(i_functions)
        # Define the voltage function
        v(t) = V_DC + quadgk(t -> -1 / C_b * i(t), 0, t)[1]

        # Initialize the search interval
        t = try
            find_zero(t -> v(t) - V_min, (interval), Bisection())
        catch
            NaN
        end

        # Check if t is NaN or Inf and replace with regression value if needed
        if isnan(t) || isinf(t)
            valid_times = filter(x -> x !== nothing && !isnan(x) && !isinf(x), t_vmin[1:idx-1])
            if !isempty(valid_times)
            t = valid_times[end] + (idx > 1 ? valid_times[end] - valid_times[end-1] : 0)  # Linear increase
            else
            t = Inf  # Default to Inf if no valid past results exist
            end
        end

        if t >= 0.9 * interval[2]
            interval = (0.0, interval[2] * 1.1)
            println(interval)
        end
        t_vmin[idx] = t
    end

    return t_vmin
end

"""
    Imax(i_functions, I_max, t_int, n)

Calculates the time at which each current function in `i_functions` first exceeds the specified maximum current `I_max`.

# Arguments
- `i_functions::Vector{Function}`: A vector of current functions, each taking a time `t` as input and returning the current at that time.
- `I_max::Float64`: The maximum current threshold to check against.
- `t_int::Float64`: The total time interval over which to evaluate the current functions.
- `n::Int`: The number of time steps to divide the interval `[0, t_int]` into.

# Returns
- `Vector{Union{Float64, Nothing}}`: A vector where each element corresponds to the time at which the respective current function first exceeds `I_max`. If the current function never exceeds `I_max`, the value will be `Inf`.

# Notes
- The time interval `[0, t_int]` is divided into `n+1` equally spaced points.
- If a current function does not exceed `I_max` at any point in the interval, the corresponding result will be `Inf`.
"""
function Imax(i_functions, I_max, t_int, n)
    time = range(0, t_int, length=n+1)
    t_imax = Vector{Union{Float64, Nothing}}(undef, length(i_functions))

    for (idx, i) in enumerate(i_functions)
        t_imax[idx] = findfirst(t -> i(t) >= I_max, time) |> x -> x !== nothing ? time[x] : Inf
    end

    return t_imax
end

"""
    I²t_time(i_functions, t_int, I²t)

Calculates the time at which the integral of the square of current functions (`i_functions`) exceeds a given threshold (`I²t`).

# Arguments
- `i_functions::Vector{Function}`: A vector of current functions, where each function represents current as a function of time.
- `t_int::Number`: The maximum time interval to consider for the integration.
- `I²t::Number`: The threshold value for the integral of the square of the current.

# Returns
- `Vector{Union{Number, Inf}}`: A vector containing the time values at which the integral of the square of each current function reaches or exceeds the threshold `I²t`. If the threshold is not reached within the given time interval, `Inf` is returned for that function.

# Notes
- The function uses numerical integration (`quadgk`) to compute the integral of the square of the current functions.
- The `find_zero` function with the `Bisection()` method is used to find the time at which the integral equals the threshold.
- The `unit` and `ustrip` functions are used to handle unit-aware calculations.

# Example
"""
function I²t_time(i_functions, t_int, I²t)
    t_I²t = []
    i²t_functions = []

    for i in i_functions
        i²t(t) = quadgk(t -> i(t)^2, 0, t)[1]
        push!(i²t_functions, i²t)
    end

    for i²t in i²t_functions
        if i²t(ustrip(t_int)) >= I²t
            t = find_zero(t -> i²t(t) - I²t, (0, ustrip(t_int)), Bisection())
            push!(t_I²t, t * unit(t_int))
        else
            push!(t_I²t, Inf)
        end
    end
    return t_I²t
end

"""
    fault_clear_prob(L_vector::StepRange, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ; return_mean::Bool=true)

Calculates the probability of fault clearance in a DC system based on various parameters.

# Arguments
- `L_vector::StepRange`: A range of inductance values.
- `L_p`: Inductance parameter.
- `C_b`: Capacitance of the DC system.
- `V_DC`: DC voltage level.
- `I_max`: Maximum allowable current.
- `V_min`: Minimum allowable voltage.
- `n`: A scaling factor for current.
- `t_max`: Maximum time for fault clearance.
- `μ`: Mean adjustment parameter for the log-normal distribution.
- `λ`: Scale parameter for the log-normal distribution.
- `return_mean::Bool=true`: If `true`, returns the mean of the cumulative distribution function (CDF); otherwise, returns the CDF values.

# Returns
- If `return_mean` is `true`, returns the mean probability of fault clearance.
- If `return_mean` is `false`, returns the CDF values for fault clearance times.

# Notes
- The function uses helper functions `sc_values`, `Vmin`, and `Imax` to compute intermediate values.
- The fault clearance times are modeled using a log-normal distribution.
"""
function fault_clear_prob(L_vector::StepRange, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ; return_mean::Bool=true)
    i_func = sc_values(ustrip(L_vector), L_p, C_b, V_DC)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    t_imax = Imax(i_func, I_max, t_vmin[end], n)
    dist = LogNormal.(log.(t_imax .+ μ) .* unit(t_max), λ * unit(t_max))

    return return_mean ? mean(_MSS.cdf.(dist, t_vmin .* unit(t_max))) : _MSS.cdf.(dist, t_vmin .* unit(t_max))
end

"""
    fault_clear_prob(L_max::Quantity, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ; return_mean::Bool=true)

Calculates the probability of clearing a fault in a DC system based on various system parameters.

# Arguments
- `L_max::Quantity`: The maximum inductance in the system.
- `L_p`: The inductance per unit length or a related parameter.
- `C_b`: The capacitance of the DC bus.
- `V_DC`: The nominal DC voltage of the system.
- `I_max`: The maximum allowable current in the system.
- `V_min`: The minimum allowable voltage in the system.
- `n`: A parameter related to the system's fault-clearing characteristics.
- `t_max`: The maximum time duration for fault clearing.
- `μ`: A parameter affecting the mean of the log-normal distribution.
- `λ`: A parameter affecting the standard deviation of the log-normal distribution.
- `return_mean::Bool=true`: If `true`, returns the mean of the cumulative distribution function (CDF); otherwise, returns the CDF value.

# Returns
- If `return_mean` is `true`, returns the mean value of the CDF for the fault-clearing time.
- If `return_mean` is `false`, returns the CDF value for the fault-clearing time.

# Notes
- The function uses `sc_values` to compute the fault current profile and `Vmin`/`Imax` to determine the critical times for voltage and current thresholds.
- A log-normal distribution is used to model the fault-clearing time, parameterized by the provided `μ` and `λ` values.
"""
function fault_clear_prob(L_max::Quantity, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ; return_mean::Bool=true)
    i_func = sc_values(ustrip(L_max), L_p, C_b, V_DC)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    t_imax = Imax(i_func, I_max, t_vmin[end], n)
    dist = LogNormal(log(t_imax[end] + μ) * unit(t_max[end]), λ * unit(t_max[end]))

    return return_mean ? mean(_MSS.cdf(dist, t_vmin[end] * unit(t_max))) : _MSS.cdf(dist, t_vmin[end] * unit(t_max))
end

"""
    fault_clear_prob_fuse(L_vector::StepRange, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ; return_mean::Bool=true)

Calculates the probability of fault clearance for a DC system protected by a fuse.

# Arguments
- `L_vector::StepRange`: A range of inductance values to evaluate.
- `L_p`: The inductance of the DC system.
- `C_b`: The capacitance of the DC system.
- `V_DC`: The DC voltage of the system.
- `I²t`: The energy let-through rating of the fuse (in terms of current squared times time).
- `V_min`: The minimum voltage threshold for the system.
- `n`: A parameter related to the system's configuration (not explicitly used in the function).
- `t_max`: The maximum time range for evaluation.
- `λ`: A scaling factor for the log-normal distribution.
- `return_mean::Bool=true`: If `true`, returns the mean probability of fault clearance; otherwise, returns the full probability distribution.

# Returns
- If `return_mean` is `true`, returns the mean probability of fault clearance as a scalar.
- If `return_mean` is `false`, returns an array of probabilities corresponding to the time values.

# Details
1. Computes the short-circuit current profile using `sc_values`.
2. Determines the time at which the voltage drops below the minimum threshold using `Vmin`.
3. Calculates the time required for the fuse to clear the fault based on the `I²t` rating using `I²t_time`.
4. Models the fault clearance time as a log-normal distribution with parameters derived from the calculated times and the scaling factor `λ`.
5. Computes the cumulative distribution function (CDF) of the log-normal distribution at the fault clearance times.

# Notes
- The function relies on helper functions `sc_values`, `Vmin`, and `I²t_time`, as well as the `_MSS.cdf` method for cumulative distribution computation.
- The `unit` function is used to ensure consistent units for time and other parameters.

"""
function fault_clear_prob_fuse(L_vector::StepRange, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ; return_mean::Bool=true)
    i_func = sc_values(ustrip(L_vector), L_p, C_b, V_DC)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    t_i²t = I²t_time(i_func, t_vmin[end], I²t)
    dist = LogNormal.(log.(t_i²t) .* unit(t_max), λ * unit(t_max))

    return return_mean ? mean(_MSS.cdf.(dist, t_vmin .* unit(t_max))) : _MSS.cdf.(dist, t_vmin .* unit(t_max))
end

"""
    fault_clear_prob_fuse(L_max::Quantity, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ; return_mean::Bool=true)

Calculates the probability of fault clearance for a DC system protected by a fuse.

# Arguments
- `L_max::Quantity`: The maximum inductance in the system.
- `L_p`: The parasitic inductance of the system.
- `C_b`: The capacitance of the DC bus.
- `V_DC`: The nominal DC voltage of the system.
- `I²t`: The energy let-through rating of the fuse (I²t value).
- `V_min`: The minimum voltage threshold for the system.
- `n`: The number of parallel branches in the system.
- `t_max`: The maximum time duration for fault clearance.
- `λ`: The scale parameter for the log-normal distribution.
- `return_mean::Bool=true`: If `true`, returns the mean of the cumulative distribution function (CDF); otherwise, returns the CDF itself.

# Returns
- If `return_mean` is `true`, returns the mean value of the fault clearance probability.
- If `return_mean` is `false`, returns the cumulative distribution function (CDF) of the fault clearance probability.

# Notes
- The function uses the `sc_values` function to compute the short-circuit current profile.
- The `Vmin` function determines the time at which the voltage drops below the minimum threshold.
- The `I²t_time` function calculates the time at which the fuse's I²t limit is reached.
- A log-normal distribution is used to model the fault clearance time.

# Dependencies
This function relies on `_MSS.cdf` for computing the cumulative distribution function of the log-normal distribution.
"""
function fault_clear_prob_fuse(L_max::Quantity, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ; return_mean::Bool=true)
    i_func = sc_values(ustrip(L_max), L_p, C_b, V_DC)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    t_i²t = I²t_time(i_func, t_vmin[end], I²t)
    dist = LogNormal(log(t_i²t[end]) * unit(t_max[end]), λ * unit(t_max[end]))

    return return_mean ? mean(_MSS.cdf(dist, t_vmin[end] * unit(t_max))) : _MSS.cdf(dist, t_vmin[end] * unit(t_max))
end