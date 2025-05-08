################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
using MultiStateSystems
using Unitful
using SpecialFunctions

"""
    Create_std(P, Pc, λᶜ)

Creates a standard state transition diagram (STD) with specified parameters.

# Arguments
- `P::Float64`: Probability parameter for certain transitions.
- `Pc::Float64`: Probability parameter for converter-related transitions.
- `λᶜ::Float64`: Rate parameter for exponential distributions.

# Returns
- `STD`: A state transition diagram object with defined states and transitions.

# Details
- The function initializes an `STD` object and adds states with specified names, power levels, and initial probabilities.
- Transitions between states are defined with various probability distributions:
  - `Exponential`: For transitions with rate `1/λᶜ` and probabilities based on `P`.
  - `LogNormal`: For transitions with log-normal distributions.
  - `Weibull`: For transitions with Weibull distributions, parameterized by `P` and `Pc`.

# Notes
- Units are specified for power (`MW`) and time (`d`, `yr`).
- Ensure that the `u"..."` macro for units is properly defined in the environment.
"""
function Create_std(P, Pc, λᶜ)    
    stdᶜᵇ = STD()
    # add the states to the std
    add_states!(stdᶜᵇ, name  = ["A", "V1", "U2","V1", "U2","V1", "U2"],
        power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW"],
        init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(2,1),(3,1),(1,4),(1,5),(4,1),(5,1),(1,6),(1,7),(6,1),(7,1)],
        distr = [   Exponential(1/(λᶜ), P/3),
                    Exponential(1/(λᶜ), (1-P)/3),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    Weibull(21.0u"yr", 3.51, P/3), # Cable
                    Weibull(21.0u"yr", 3.51, (1-P)/3), # Cable
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    Weibull(15.0u"yr", 3.81, Pc/3), # Converter
                    Weibull(15.0u"yr", 3.81, (1-Pc)/3), # Converter
                    LogNormal(log(10.0)u"d", 0.2u"d"),
                    LogNormal(log(10.0)u"d", 0.2u"d")])
    return stdᶜᵇ
end

"""
    Create_std(P, λᶜ)

Creates a standard (STD) object with predefined states and transitions.

# Arguments
- `P::Float64`: A parameter that determines the probability distribution of certain transitions.
- `λᶜ::Float64`: A parameter used to define the rate of exponential distributions.

# Returns
- `STD`: A standard object with states and transitions added.

# Details
- The function initializes an `STD` object and adds states with specific names, power levels, and initial probabilities.
- Transitions between states are defined with various probability distributions:
  - `Exponential`: Used for transitions with rates dependent on `λᶜ` and `P`.
  - `LogNormal`: Used for transitions with a log-normal distribution.
  - `Weibull`: Used for transitions with a Weibull distribution, representing cable-related transitions.

# Notes
- The `power` parameter for states includes units (e.g., `u"MW"` for megawatts).
- The `distr` parameter for transitions includes units (e.g., `u"d"` for days, `u"yr"` for years).
"""
function Create_std(P, λᶜ)    
    stdᶜᵇ = STD()
    # add the states to the std
    add_states!(stdᶜᵇ, name  = ["A", "V1", "U2","V1", "U2"],
        power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW"],
        init  = [1.0, 0.0, 0.0, 0.0, 0.0])
    add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(2,1),(3,1),(1,4),(1,5),(4,1),(5,1)],
        distr = [   Exponential(1/(λᶜ), P/2),
                    Exponential(1/(λᶜ), (1-P)/2),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    Weibull(21.0u"yr", 3.51, P/2), # Cable
                    Weibull(21.0u"yr", 3.51, (1-P)/2), # Cable
                    LogNormal(log(14.0)u"d", 0.1u"d"),
                    LogNormal(log(14.0)u"d", 0.1u"d")])
    return stdᶜᵇ
end

"""
    Create_std(P, Pc, λᶜ; type = "SemiMarkov", P_zone = 1)

Creates a standard state transition diagram (STD) based on the specified parameters.

# Arguments
- `P::Float64`: The probability parameter for the system.
- `Pc::Float64`: The probability parameter specific to the converter.
- `λᶜ::Float64`: The rate parameter for exponential distributions.
- `type::String`: The type of STD to create. Can be `"SemiMarkov"` or `"Markov"`. Defaults to `"SemiMarkov"`.
- `P_zone::Float64`: The scaling factor for the zone. Defaults to `1`.

# Returns
- `stdᶜᵇ`: The created STD object.

# Details
- For `type = "SemiMarkov"`, the function defines states and transitions using a mix of `Exponential`, `LogNormal`, and `Weibull` distributions.
- For `type = "Markov"`, the function defines states and transitions using only `Exponential` distributions.
- The states represent different operational and failure modes of the system, while the transitions define the probabilities and rates of moving between these states.

# Errors
- Throws an error if an invalid `type` is provided.
"""
function Create_std(P, Pc, λᶜ; type = "SemiMarkov", P_zone = 1)    
    stdᶜᵇ = STD()    
    if type == "SemiMarkov"
    # add the states to the std
    add_states!(stdᶜᵇ, name  = ["A", "V1", "U2", "U3", "V1", "U2", "U3", "V1", "U2", "U3"],
        power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW"],
        init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    
        add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(1,4),(2,1),(3,1),(4,1),(1,5),(1,6),(1,7),(5,1),(6,1),(7,1),(1,8),(1,9),(1,10),(8,1),(9,1),(10,1)],
            distr = [   Exponential(1/(λᶜ), (P*P_zone)/3),
                        Exponential(1/(λᶜ), (1-P)*P_zone/3),
                        Exponential(1/(λᶜ), (1-P_zone)/3),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        Weibull(21.0u"yr", 3.51, P*P_zone/3), # Cable
                        Weibull(21.0u"yr", 3.51, (1-P)*P_zone/3), # Cable
                        Weibull(21.0u"yr", 3.51, (1-P_zone)/3),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        Weibull(15.0u"yr", 3.81, Pc*P_zone/3), # Converter
                        Weibull(15.0u"yr", 3.81, (1-Pc)*P_zone/3), # Converter
                        Weibull(15.0u"yr", 3.81, (1-P_zone)/3),
                        LogNormal(log(10.0)u"d", 0.2u"d"),
                        LogNormal(log(10.0)u"d", 0.2u"d"),
                        LogNormal(log(10.0)u"d", 0.2u"d")])
    elseif type == "Markov"
        add_states!(stdᶜᵇ, name  = ["A", "V1", "U2", "V2", "U3", "V3", "V1", "U2", "V2", "U3", "V3", "V1", "U2", "V2", "U3", "V3"],
        power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW"],
        init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  
    add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(1,5),(3,4),(5,6),(2,1),(4,1),(6,1),(1,7),(1,8),(1,10),(8,9),(10,11),(7,1),(9,1),(11,1),(1,12),(1,13),(1,15),(13,14),(15,16),(12,1),(14,1),(16,1)], 
            distr = [   Exponential(1/(λᶜ), (P*P_zone)/3),
                        Exponential(1/(λᶜ), (1-P)*P_zone/3),
                        Exponential(1/(λᶜ), (1-P_zone)/3),
                        Exponential(10.0u"hr"),
                        Exponential(10.0u"hr"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(21.0u"yr", P*P_zone/3), # Cable
                        Exponential(21.0u"yr", (1-P)*P_zone/3), # Cable
                        Exponential(21.0u"yr", (1-P_zone)/3),
                        Exponential(10.0u"hr"),
                        Exponential(10.0u"hr"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(15.0u"yr", Pc*P_zone/3), # Converter
                        Exponential(15.0u"yr", (1-Pc)*P_zone/3), # Converter
                        Exponential(15.0u"yr", (1-P_zone)/3),
                        Exponential(10.0u"hr"),
                        Exponential(10.0u"hr"),
                        Exponential(10.0u"d"),
                        Exponential(10.0u"d"),
                        Exponential(10.0u"d")])
    else
        error("Invalid type. Use 'SemiMarkov' or 'Markov'.")
    end
    return stdᶜᵇ
end

"""
    Create_std(P, λᶜ; type = "SemiMarkov", P_zone = 1)

Creates a standard state transition diagram (STD) based on the specified parameters.

# Arguments
- `P::Float64`: A parameter representing the probability of a specific event.
- `λᶜ::Float64`: A parameter representing the rate of occurrence for certain transitions.
- `type::String`: (Optional) Specifies the type of STD to create. Can be either `"SemiMarkov"` or `"Markov"`. Defaults to `"SemiMarkov"`.
- `P_zone::Float64`: (Optional) A scaling factor for the probabilities. Defaults to `1`.

# Returns
- `STD`: A state transition diagram object with states and transitions defined based on the input parameters.

# Details
- If `type` is `"SemiMarkov"`, the function creates a semi-Markov STD with specific states and transitions, using distributions such as `Exponential`, `LogNormal`, and `Weibull`.
- If `type` is `"Markov"`, the function creates a Markov STD with a different set of states and transitions, using `Exponential` distributions.

# Errors
- Throws an error if `type` is not `"SemiMarkov"` or `"Markov"`.
"""
function Create_std(P, λᶜ; type = "SemiMarkov", P_zone = 1)    
    stdᶜᵇ = STD()    
    if type == "SemiMarkov"
        # add the states to the std
        add_states!(stdᶜᵇ, name  = ["A", "V1", "U2", "U3", "V1", "U2", "U3"],
            power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW"],
            init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])    
        add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(1,4),(2,1),(3,1),(4,1),(1,5),(1,6),(1,7),(5,1),(6,1),(7,1)],
            distr = [   Exponential(1/(λᶜ), (P*P_zone)/2),
                        Exponential(1/(λᶜ), (1-P)*P_zone/2),
                        Exponential(1/(λᶜ), (1-P_zone)/2),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        Weibull(21.0u"yr", 3.51, P*P_zone/2), # Cable
                        Weibull(21.0u"yr", 3.51, (1-P)*P_zone/2), # Cable
                        Weibull(21.0u"yr", 3.51, (1-P_zone)/2),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),
                        LogNormal(log(14.0)u"d", 0.1u"d"),])
    elseif type == "Markov"
        add_states!(stdᶜᵇ, name  = ["A", "V1", "U2", "V2", "U3", "V3", "V1", "U2", "V2", "U3", "V3"],
            power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW"],
            init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  
        add_transitions!(stdᶜᵇ, states = [(1,2),(1,3),(1,5),(3,4),(5,6),(2,1),(4,1),(6,1),(1,7),(1,8),(1,10),(8,9),(10,11),(7,1),(9,1),(11,1)], 
            distr = [   Exponential(1/(λᶜ), (P*P_zone)/2),
                        Exponential(1/(λᶜ), (1-P)*P_zone/2),
                        Exponential(1/(λᶜ), (1-P_zone)/2),
                        Exponential(10.0u"hr"),
                        Exponential(10.0u"hr"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(21.0u"yr", P*P_zone/2), # Cable
                        Exponential(21.0u"yr", (1-P)*P_zone/2), # Cable
                        Exponential(21.0u"yr", (1-P_zone)/2),
                        Exponential(10.0u"hr"),
                        Exponential(10.0u"hr"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d"),
                        Exponential(14.0u"d")])
    else 
        error("Invalid type. Use 'SemiMarkov' or 'Markov'.")
    end
    return stdᶜᵇ
end

"""
    fill_std(L_tot, P, λᶜ; type, Pc=nothing, P_zone=1)

Creates a nested dictionary (`std`) containing standardized data for each key in `L_tot`.

# Arguments
- `L_tot::Dict`: A dictionary where each key maps to another dictionary (`L_c`) containing data to be processed.
- `P::Dict`: A dictionary containing parameters associated with the keys in `L_tot`.
- `λᶜ::Number`: A scaling factor applied to the maximum value of each entry in `L_c`.
- `type`: A keyword argument specifying the type of standardization to be applied.
- `Pc::Union{Dict, Nothing}` (optional): An optional dictionary containing additional parameters for standardization. Defaults to `nothing`.
- `P_zone::Union{Number, Nothing}` (optional): An optional parameter specifying the zone for standardization. Defaults to `1`.

# Returns
- `std::Dict`: A nested dictionary where each key corresponds to a key in `L_tot`, and each value is a dictionary containing standardized data.

# Details
- For each key in `L_tot`, the function iterates over its associated dictionary (`L_c`).
- Depending on the presence of `Pc` and `P_zone`, the function calls `Create_std` with different combinations of arguments.
- The `Create_std` function is expected to handle the actual standardization process.

# Notes
- The `Create_std` function must be defined elsewhere in the codebase.
- The `type` and `P_zone` arguments are passed to `Create_std` when provided.
"""
function fill_std(L_tot, P, λᶜ; type, Pc = nothing, P_zone = 1)
    std = Dict()
    for (key, L_c) in L_tot
        std_i = Dict()
        for (key_a, value_a) in L_c
            if isnothing(P_zone)
                if isnothing(Pc)
                    std_i[key_a] = Create_std(P[key][key_a], λᶜ * maximum(value_a))
                else
                    std_i[key_a] = Create_std(P[key][key_a], Pc[key][key_a], λᶜ * maximum(value_a))
                end
            else
                if isnothing(Pc)
                    std_i[key_a] = Create_std(P[key][key_a], λᶜ * maximum(value_a); type = type, P_zone = P_zone)
                else
                    std_i[key_a] = Create_std(P[key][key_a], Pc[key][key_a], λᶜ * maximum(value_a); type = type, P_zone = P_zone)
                end
            end
        end
        std[key] = std_i
    end
    return std
end

"""
    solve_std_dict!(std::AbstractDict, cls, tsim::Number, dt::Number)

Recursively solves a nested dictionary of standard objects (`std`) using the provided class (`cls`), simulation time (`tsim`), 
and time step (`dt`). 

# Arguments
- `std::AbstractDict`: A dictionary containing standard objects or nested dictionaries of standard objects.
- `cls`: The class or configuration used for solving the standard objects.
- `tsim::Number`: The total simulation time.
- `dt::Number`: The time step for the simulation.

# Behavior
- If a value in the dictionary is itself a dictionary, the function calls itself recursively to process the nested dictionary.
- If a value is not a dictionary, it is assumed to be a standard object and is solved using `_MSS.solve!`.

# Notes
- The `_MSS.solve!` function is expected to handle the actual solving process for the standard objects.
- This function modifies the input dictionary `std` in place.

"""
function solve_std_dict!(std, cls, tsim, dt)
    for (key, std_val) in std
        if isa(std_val, AbstractDict)
            solve_std_dict!(std_val, cls, tsim, dt)
        else
            # Do something with the value
            _MSS.solve!(std_val, cls; tsim, dt)
            # println("Solved STD for $key.")
        end
    end
end