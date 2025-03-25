################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
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
    sc_values(L_vector::Union{AbstractVector, Int}, L_p::Float64, C_b::Float64, V_DC::Int64)

Calculate short-circuit current functions for a given set of cable lengths.

# Arguments
- `L_vector::Union{AbstractVector, Int}`: A vector of cable lengths or a single integer representing the cable length(s).
- `L_p::Float64`: The inductance of the power source in Henrys (H).
- `C_b::Float64`: The capacitance of the cable in Farads (F).
- `V_DC::Int64`: The DC voltage in Volts (V).

# Returns
- `Vector{Function}`: A vector of functions, each representing the short-circuit current as a function of time for the corresponding cable length in `L_vector`.

# Details
The function calculates the short-circuit current for each cable length in `L_vector` based on the given parameters. It handles both underdamped and overdamped cases:
- **Underdamped**: When the resistance squared is less than four times the inductance divided by the capacitance.
- **Overdamped**: When the resistance squared is greater than or equal to four times the inductance divided by the capacitance.

The resulting current functions are stored in a vector and returned.
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
            i = t -> V_DC / (ω * L_total) * exp(-R / (2 * L_total) * t) * sinh(ω * t)
        end
        i_functions[idx] = i
    end
    return i_functions
end

"""
    Imax(i_functions, I_max, t_int, n)

Calculate the time instances at which the current exceeds a specified maximum value.

# Arguments
- `i_functions::Vector{Function}`: A vector of current functions, each representing the current as a function of time.
- `I_max::Number`: The maximum allowable current value.
- `t_int::Number`: The total time interval to consider.
- `n::Int`: The number of time steps to divide the interval into.

# Returns
- `t_imax::Vector{Number}`: A vector of time instances at which the current exceeds `I_max` for each function in `i_functions`.

"""
function Imax(i_functions, I_max, t_int, n)
    time = range(0, t_int, length = n+1)
    t_imax = []
    # I max criterium
    for i in i_functions
        for t in time
            i(t) >= I_max ? (push!(t_imax, t), break) : ~ ;            
        end
    end
    return t_imax
end

function Vmin(i_functions, C_b, t_int, V_DC, V_min)
    v_functions = []
    t_vmin = []
    for i in i_functions
        v(t) = V_DC + quadgk(t -> -1/C_b * i(t), 0, t)[1]
        push!(v_functions, v)
    end

    # V min criterium
    for v in v_functions
        t = find_zero(t -> v(t) - V_min, (0, ustrip(t_int)), Bisection())
        if t <= ustrip(t_int)
            push!(t_vmin, t)
        else
            push!(t_vmin, nothing)
        end
    end

    return t_vmin
end

"""
    I²t_time(i_functions, t_int, I²t)

Calculate the time at which the integral of the square of the current (`I²t`) reaches a specified value for a set of current functions.

# Arguments
- `i_functions::Vector{Function}`: A vector of functions representing the current as a function of time.
- `t_int::Number`: The time interval over which to evaluate the integral.
- `I²t::Number`: The target value of the integral of the square of the current.

# Returns
- `t_I²t::Vector{Union{Number, Inf}}`: A vector of times at which the integral of the square of the current reaches the specified value. If the integral does not reach the specified value within the given time interval, `Inf` is returned for that function.

"""
function I²t_time(i_functions, t_int, I²t)
    t_I²t  = []
    i²t_functions = []

    for i in i_functions
        i²t(t) = quadgk((t -> i(t)^2), 0, t)[1]
        push!(i²t_functions, i²t)
    end

    # I²t criterium
    for i²t in i²t_functions
        a, b = 0, ustrip(t_int)
        if i²t(a) * i²t(b) < 0
            t = find_zero(t -> i²t(t) - I²t, (a, b), Bisection())
            t_with_units = t * ustrip(unit(t_int))
            push!(t_I²t, t_with_units)
        else
            push!(t_I²t, Inf)
        end
    end
    
    if length(t_I²t) == 0
        t_I²t = [Inf]
    elseif length(t_I²t) < length(i_functions)
        for i in 1:(length(i_functions)-length(t_I²t))-1
        push!(t_I²t, Inf)
        end
    end

    return t_I²t
end


"""
    fault_clear_prob(L_vector::StepRange, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)

Calculate the probability of fault clearance for a given range of inductance values.

# Arguments
- `L_vector::StepRange`: Range of inductance values.
- `L_p`: Inductance parameter.
- `C_b`: Capacitance value.
- `V_DC`: DC voltage.
- `I_max`: Maximum current.
- `V_min`: Minimum voltage.
- `n`: Number of samples.
- `t_max`: Maximum time.
- `μ`: Mean of the log-normal distribution.
- `λ`: Standard deviation of the log-normal distribution.

# Returns
- `Float64`: Mean probability of fault clearance.
"""
function fault_clear_prob(L_vector::StepRange, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
    i_func = sc_values(ustrip(L_vector), L_p, C_b, V_DC)
    t_imax = Imax(i_func, I_max, ustrip(t_max), n)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    dist = LogNormal.(log.(t_imax .+ μ) .* unit(t_max), λ * unit(t_max))
    return mean(_MSS.cdf.(dist, t_vmin .* unit(t_max)))
end

"""
    fault_clear_prob(L_max::Quantity, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)

Calculate the probability of fault clearance for a given maximum inductance value.

# Arguments
- `L_max::Quantity`: Maximum inductance value.
- `L_p`: Inductance parameter.
- `C_b`: Capacitance value.
- `V_DC`: DC voltage.
- `I_max`: Maximum current.
- `V_min`: Minimum voltage.
- `n`: Number of samples.
- `t_max`: Maximum time.
- `μ`: Mean of the log-normal distribution.
- `λ`: Standard deviation of the log-normal distribution.

# Returns
- `Float64`: Mean probability of fault clearance.
"""
function fault_clear_prob(L_max::Quantity, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
    i_func = sc_values(ustrip(L_max), L_p, C_b, V_DC)
    t_imax = Imax(i_func, I_max, ustrip(t_max), n)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    dist = LogNormal(log(t_imax[end] + μ) * unit(t_max[end]), λ * unit(t_max[end]))
    return mean(_MSS.cdf(dist, t_vmin[end] * unit(t_max)))
end

"""
    fault_clear_prob_fuse(L_vector::StepRange, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ)

Calculate the probability of fault clearance for a fuse given a range of inductance values.

# Arguments
- `L_vector::StepRange`: Range of inductance values.
- `L_p`: Inductance parameter.
- `C_b`: Capacitance value.
- `V_DC`: DC voltage.
- `I²t`: Energy let-through of the fuse.
- `V_min`: Minimum voltage.
- `n`: Number of samples.
- `t_max`: Maximum time.
- `λ`: Standard deviation of the log-normal distribution.

# Returns
- `Float64`: Mean probability of fault clearance.
"""
function fault_clear_prob_fuse(L_vector::StepRange, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ)
    i_func = sc_values(ustrip(L_vector), L_p, C_b, V_DC)
    t_i²t = I²t_time(i_func, ustrip(t_max), I²t)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    dist = LogNormal.(log.(t_i²t) .* unit(t_max), λ * unit(t_max))
    return mean(_MSS.cdf.(dist, t_vmin .* unit(t_max)))
end

"""
    fault_clear_prob_fuse(L_max::Quantity, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ)

Calculate the probability of fault clearance for a fuse given a maximum inductance value.

# Arguments
- `L_max::Quantity`: Maximum inductance value.
- `L_p`: Inductance parameter.
- `C_b`: Capacitance value.
- `V_DC`: DC voltage.
- `I²t`: Energy let-through of the fuse.
- `V_min`: Minimum voltage.
- `n`: Number of samples.
- `t_max`: Maximum time.
- `λ`: Standard deviation of the log-normal distribution.

# Returns
- `Float64`: Mean probability of fault clearance.
"""
function fault_clear_prob_fuse(L_max::Quantity, L_p, C_b, V_DC, I²t, V_min, n, t_max, λ)
    i_func = sc_values(ustrip(L_max), L_p, C_b, V_DC)
    t_i²t = I²t_time(i_func, ustrip(t_max), I²t)
    t_vmin = Vmin(i_func, C_b, ustrip(t_max), V_DC, V_min)
    dist = LogNormal(log(t_i²t[end]) * unit(t_max[end]), λ * unit(t_max[end]))
    return mean(_MSS.cdf(dist, t_vmin[end] * unit(t_max)))
end