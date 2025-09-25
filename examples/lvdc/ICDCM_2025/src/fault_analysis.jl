################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Fault Analysis Module

This module contains functions for DC fault current analysis and protection device modeling.
It provides the core fault analysis capabilities used in the LVDC short-circuit examples.

## Main Functions
- `calculate_short_circuit_current`: Compute fault current profiles
- `calculate_voltage_drop_time`: Find time when voltage drops below threshold
- `calculate_overcurrent_time`: Find time when current exceeds threshold  
- `calculate_energy_time`: Find time when I²t energy limit is reached
- `calculate_fault_clearing_probability`: Compute protection device clearing probabilities
"""

using MultiStateSystems
using Unitful
using Statistics
using QuadGK
using Roots

# Default cable parameters
const CABLE_INDUCTANCE_PER_METER = 2.94e-7  # H/m
const CABLE_RESISTANCE_PER_METER = 3.08e-3  # Ω/m

"""
    calculate_short_circuit_current(cable_lengths, L_p::Float64, C_b::Float64, V_DC::Int)

Calculate short-circuit current functions for given cable lengths.

# Arguments
- `cable_lengths`: Vector of cable lengths (m) or single length
- `L_p::Float64`: Parasitic inductance (H)
- `C_b::Float64`: Bus capacitance (F)
- `V_DC::Int`: DC voltage (V)

# Returns
- `Vector{Function}`: Current functions i(t) for each cable length
"""
function calculate_short_circuit_current(cable_lengths, L_p::Float64, C_b::Float64, V_DC::Int)
    # Handle single length input
    lengths = isa(cable_lengths, Number) ? [cable_lengths] : cable_lengths
    
    current_functions = Vector{Function}(undef, length(lengths))
    
    for (idx, L) in enumerate(lengths)
        R = 2 * CABLE_RESISTANCE_PER_METER * L
        L_total = 2 * CABLE_INDUCTANCE_PER_METER * L + L_p
        
        if R^2 < 4 * L_total / C_b  # Underdamped response
            ω = sqrt(1 / (L_total * C_b) - (R / (2 * L_total))^2)
            current_functions[idx] = t -> V_DC / (ω * L_total) * exp(-R / (2 * L_total) * t) * sin(ω * t)
        else  # Overdamped response
            ω = sqrt((R / (2 * L_total))^2 - 1 / (L_total * C_b))
            current_functions[idx] = t -> begin
                exp_term = exp(-R / (2 * L_total) * t)
                # Prevent numerical overflow
                if ω * t > 700
                    sinh_term = 2 * exp(ω * t - 700) * exp(700)
                else
                    sinh_term = exp(ω * t) - exp(-ω * t)
                end
                V_DC / (2 * ω * L_total) * exp_term * sinh_term
            end
        end
    end
    
    return current_functions
end

"""
    calculate_voltage_drop_time(current_functions, C_b::Float64, t_max, V_DC::Int, V_min::Float64)

Calculate the time when DC bus voltage drops below minimum threshold.

# Arguments
- `current_functions`: Vector of current functions
- `C_b::Float64`: Bus capacitance (F)
- `t_max`: Maximum time to consider
- `V_DC::Int`: Initial DC voltage (V)
- `V_min::Float64`: Minimum voltage threshold (V)

# Returns
- `Vector`: Times when voltage drops below V_min for each current function
"""
function calculate_voltage_drop_time(current_functions, C_b::Float64, t_max, V_DC::Int, V_min::Float64)
    t_vmin = Vector{Union{Float64, Nothing}}(undef, length(current_functions))
    interval = (0.0, ustrip(t_max))
    
    for (idx, i_func) in enumerate(current_functions)
        # Voltage function: v(t) = V_DC - ∫(i(t)/C_b)dt from 0 to t
        voltage_function(t) = V_DC + quadgk(τ -> -1/C_b * i_func(τ), 0, t)[1]
        
        try
            t = find_zero(t -> voltage_function(t) - V_min, interval, Bisection())
            t_vmin[idx] = t
        catch
            # Handle cases where no solution is found
            valid_times = filter(x -> x !== nothing && !isnan(x) && !isinf(x), t_vmin[1:idx-1])
            if !isempty(valid_times)
                t_vmin[idx] = valid_times[end] + (idx > 1 ? valid_times[end] - valid_times[end-1] : 0)
            else
                t_vmin[idx] = Inf
            end
        end
        
        # Expand search interval if needed
        if !isnothing(t_vmin[idx]) && t_vmin[idx] >= 0.9 * interval[2]
            interval = (0.0, interval[2] * 1.1)
        end
    end
    
    return t_vmin
end

"""
    calculate_overcurrent_time(current_functions, I_max::Int, t_max, n::Int)

Calculate when current functions exceed maximum current threshold.

# Arguments
- `current_functions`: Vector of current functions
- `I_max::Int`: Maximum current threshold (A)
- `t_max`: Maximum time to evaluate
- `n::Int`: Number of time points

# Returns
- `Vector`: Times when current exceeds I_max for each function
"""
function calculate_overcurrent_time(current_functions, I_max::Int, t_max, n::Int)
    time_points = range(0, ustrip(t_max), length=n+1)
    t_imax = Vector{Union{Float64, Nothing}}(undef, length(current_functions))
    
    for (idx, i_func) in enumerate(current_functions)
        t_imax[idx] = findfirst(t -> i_func(t) >= I_max, time_points) |> 
                     x -> x !== nothing ? time_points[x] : Inf
    end
    
    return t_imax
end

"""
    calculate_energy_time(current_functions, t_max, I²t_limit)

Calculate when I²t energy integral reaches the specified limit.

# Arguments
- `current_functions`: Vector of current functions
- `t_max`: Maximum time to consider
- `I²t_limit`: Energy limit (A²s)

# Returns
- `Vector`: Times when I²t limit is reached for each function
"""
function calculate_energy_time(current_functions, t_max, I²t_limit)
    t_I²t = []
    
    for i_func in current_functions
        # Energy function: E(t) = ∫(i(t)²)dt from 0 to t
        energy_function(t) = quadgk(τ -> i_func(τ)^2, 0, t)[1]
        
        if energy_function(ustrip(t_max)) >= I²t_limit
            t = find_zero(t -> energy_function(t) - I²t_limit, (0, ustrip(t_max)), Bisection())
            push!(t_I²t, t * unit(t_max))
        else
            push!(t_I²t, Inf)
        end
    end
    
    return t_I²t
end

"""
    calculate_fault_clearing_probability(cable_lengths, system_params, int_params, device_params; return_mean::Bool=true)

Calculate the probability of successful fault clearing for circuit breakers.

# Arguments
- `cable_lengths`: Cable length range or single value
- `system_params`: System parameters structure
- `int_params`: Integration parameters structure
- `device_params`: Protection device parameters (μ, λ)
- `return_mean::Bool`: Return mean probability vs full CDF

# Returns
- Fault clearing probability (mean or CDF array)
"""
function calculate_fault_clearing_probability(cable_lengths, system_params, int_params, device_params; return_mean::Bool=true)
    # Calculate short-circuit current profiles
    i_functions = calculate_short_circuit_current(
        ustrip(cable_lengths), system_params.L_p, system_params.C_b, system_params.V_DC
    )
    
    # Calculate critical times
    t_vmin = calculate_voltage_drop_time(
        i_functions, system_params.C_b, int_params.t_max, 
        system_params.V_DC, system_params.V_min
    )
    
    t_imax = calculate_overcurrent_time(
        i_functions, system_params.I_max, t_vmin[end], int_params.n
    )
    
    # Model clearing time with log-normal distribution
    μ, λ = device_params
    if isa(cable_lengths, AbstractRange) || (isa(cable_lengths, AbstractArray) && length(cable_lengths) > 1)
        dist = LogNormal.(log.(t_imax .+ μ) .* unit(int_params.t_max), λ * unit(int_params.t_max))
        cdf_values = _MSS.cdf.(dist, t_vmin .* unit(int_params.t_max))
    else
        # Single value case
        t_val = isa(t_imax, AbstractArray) ? t_imax[1] : t_imax
        t_vmin_val = isa(t_vmin, AbstractArray) ? t_vmin[1] : t_vmin
        dist = LogNormal(log(t_val + μ) * unit(int_params.t_max), λ * unit(int_params.t_max))
        cdf_values = _MSS.cdf(dist, t_vmin_val * unit(int_params.t_max))
    end
    
    return return_mean ? mean(cdf_values) : cdf_values
end

"""
    calculate_fuse_clearing_probability(cable_lengths, system_params, int_params, λ; return_mean::Bool=true)

Calculate the probability of successful fault clearing for fuses.

# Arguments
- `cable_lengths`: Cable length range or single value
- `system_params`: System parameters structure
- `int_params`: Integration parameters structure
- `λ`: Fuse parameter
- `return_mean::Bool`: Return mean probability vs full CDF

# Returns
- Fault clearing probability (mean or CDF array)
"""
function calculate_fuse_clearing_probability(cable_lengths, system_params, int_params, λ; return_mean::Bool=true)
    # Calculate short-circuit current profiles
    i_functions = calculate_short_circuit_current(
        ustrip(cable_lengths), system_params.L_p, system_params.C_b, system_params.V_DC
    )
    
    # Calculate critical times
    t_vmin = calculate_voltage_drop_time(
        i_functions, system_params.C_b, int_params.t_max,
        system_params.V_DC, system_params.V_min
    )
    
    t_i²t = calculate_energy_time(i_functions, t_vmin[end], system_params.I²t)
    
    # Model clearing time with log-normal distribution
    if isa(cable_lengths, AbstractRange) || (isa(cable_lengths, AbstractArray) && length(cable_lengths) > 1)
        dist = LogNormal.(log.(t_i²t) .* unit(int_params.t_max), λ * unit(int_params.t_max))
        cdf_values = _MSS.cdf.(dist, t_vmin .* unit(int_params.t_max))
    else
        # Single value case  
        t_val = isa(t_i²t, AbstractArray) ? t_i²t[1] : t_i²t
        t_vmin_val = isa(t_vmin, AbstractArray) ? t_vmin[1] : t_vmin
        dist = LogNormal(log(t_val) * unit(int_params.t_max), λ * unit(int_params.t_max))
        cdf_values = _MSS.cdf(dist, t_vmin_val * unit(int_params.t_max))
    end
    
    return return_mean ? mean(cdf_values) : cdf_values
end
