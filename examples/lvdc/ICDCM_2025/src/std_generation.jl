################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# State Transition Diagram Generation Module

This module provides functions for creating and solving state transition diagrams (STDs)
for LVDC system components, including cables, converters, and protection devices.

## Main Functions
- `create_component_std`: Create STD for individual components
- `generate_component_stds`: Generate STDs for all system components
- `solve_component_stds!`: Solve all generated STDs
- `create_reduced_stds`: Create reduced state representations
"""

using MultiStateSystems
using Unitful
using SpecialFunctions

"""
    create_component_std(P::Float64, Pc::Float64, λᶜ::Float64; process_type::String="SemiMarkov", P_zone::Float64=1.0)

Create a state transition diagram for a component with converter.

# Arguments
- `P::Float64`: Component reliability probability
- `Pc::Float64`: Converter reliability probability  
- `λᶜ`: Failure rate parameter
- `process_type::String`: "SemiMarkov" or "Markov"
- `P_zone::Float64`: Zone protection factor

# Returns
- `STD`: State transition diagram object
"""
function create_component_std(P::Float64, Pc::Float64, λᶜ; process_type::String="SemiMarkov", P_zone::Float64=1.0)
    std = STD()
    
    if process_type == "SemiMarkov"
        # States: A (available), V1 (voltage fault), U2/U3 (unavailable), etc.
        add_states!(std, 
            name = ["A", "V1", "U2", "U3", "V1", "U2", "U3", "V1", "U2", "U3"],
            power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW"],
            init = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        
        # Transitions with various distributions
        add_transitions!(std, 
            states = [(1,2),(1,3),(1,4),(2,1),(3,1),(4,1),(1,5),(1,6),(1,7),(5,1),(6,1),(7,1),(1,8),(1,9),(1,10),(8,1),(9,1),(10,1)],
            distr = [
                # Short-circuit transitions (exponential)
                Exponential(1/λᶜ, (P*P_zone)/3),
                Exponential(1/λᶜ, (1-P)*P_zone/3),
                Exponential(1/λᶜ, (1-P_zone)/3),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                # Cable aging transitions (Weibull)
                Weibull(21.0u"yr", 3.51, P*P_zone/3),
                Weibull(21.0u"yr", 3.51, (1-P)*P_zone/3),
                Weibull(21.0u"yr", 3.51, (1-P_zone)/3),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                # Converter transitions (Weibull)
                Weibull(15.0u"yr", 3.81, Pc*P_zone/3),
                Weibull(15.0u"yr", 3.81, (1-Pc)*P_zone/3),
                Weibull(15.0u"yr", 3.81, (1-P_zone)/3),
                LogNormal(log(10.0)u"d", 0.2u"d"),
                LogNormal(log(10.0)u"d", 0.2u"d"),
                LogNormal(log(10.0)u"d", 0.2u"d")
            ]
        )
        
    elseif process_type == "Markov"
        # Markov process with exponential transitions only
        add_states!(std,
            name = ["A", "V1", "U2", "V2", "U3", "V3", "V1", "U2", "V2", "U3", "V3", "V1", "U2", "V2", "U3", "V3"],
            power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW"],
            init = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        
        add_transitions!(std,
            states = [(1,2),(1,3),(1,5),(3,4),(5,6),(2,1),(4,1),(6,1),(1,7),(1,8),(1,10),(8,9),(10,11),(7,1),(9,1),(11,1),(1,12),(1,13),(1,15),(13,14),(15,16),(12,1),(14,1),(16,1)],
            distr = [
                Exponential(1/λᶜ, (P*P_zone)/3),
                Exponential(1/λᶜ, (1-P)*P_zone/3),
                Exponential(1/λᶜ, (1-P_zone)/3),
                Exponential(10.0u"hr"),
                Exponential(10.0u"hr"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(21.0u"yr", P*P_zone/3),
                Exponential(21.0u"yr", (1-P)*P_zone/3),
                Exponential(21.0u"yr", (1-P_zone)/3),
                Exponential(10.0u"hr"),
                Exponential(10.0u"hr"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(15.0u"yr", Pc*P_zone/3),
                Exponential(15.0u"yr", (1-Pc)*P_zone/3),
                Exponential(15.0u"yr", (1-P_zone)/3),
                Exponential(10.0u"hr"),
                Exponential(10.0u"hr"),
                Exponential(10.0u"d"),
                Exponential(10.0u"d"),
                Exponential(10.0u"d")
            ]
        )
    else
        error("Invalid process type: $process_type. Use 'SemiMarkov' or 'Markov'.")
    end
    
    return std
end

"""
    create_component_std(P::Float64, λᶜ::Float64; process_type::String="SemiMarkov", P_zone::Float64=1.0)

Create a state transition diagram for a component without converter.

# Arguments
- `P::Float64`: Component reliability probability
- `λᶜ`: Failure rate parameter
- `process_type::String`: "SemiMarkov" or "Markov"
- `P_zone::Float64`: Zone protection factor

# Returns
- `STD`: State transition diagram object
"""
function create_component_std(P::Float64, λᶜ; process_type::String="SemiMarkov", P_zone::Float64=1.0)
    std = STD()
    
    if process_type == "SemiMarkov"
        add_states!(std,
            name = ["A", "V1", "U2", "U3", "V1", "U2", "U3"],
            power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW"],
            init = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        
        add_transitions!(std,
            states = [(1,2),(1,3),(1,4),(2,1),(3,1),(4,1),(1,5),(1,6),(1,7),(5,1),(6,1),(7,1)],
            distr = [
                Exponential(1/λᶜ, (P*P_zone)/2),
                Exponential(1/λᶜ, (1-P)*P_zone/2),
                Exponential(1/λᶜ, (1-P_zone)/2),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                Weibull(21.0u"yr", 3.51, P*P_zone/2),
                Weibull(21.0u"yr", 3.51, (1-P)*P_zone/2),
                Weibull(21.0u"yr", 3.51, (1-P_zone)/2),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d"),
                LogNormal(log(14.0)u"d", 0.1u"d")
            ]
        )
        
    elseif process_type == "Markov"
        add_states!(std,
            name = ["A", "V1", "U2", "V2", "U3", "V3", "V1", "U2", "V2", "U3", "V3"],
            power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW", (Inf)u"MW", 0.0u"MW"],
            init = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        )
        
        add_transitions!(std,
            states = [(1,2),(1,3),(1,5),(3,4),(5,6),(2,1),(4,1),(6,1),(1,7),(1,8),(1,10),(8,9),(10,11),(7,1),(9,1),(11,1)],
            distr = [
                Exponential(1/λᶜ, (P*P_zone)/2),
                Exponential(1/λᶜ, (1-P)*P_zone/2),
                Exponential(1/λᶜ, (1-P_zone)/2),
                Exponential(10.0u"hr"),
                Exponential(10.0u"hr"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(21.0u"yr", P*P_zone/2),
                Exponential(21.0u"yr", (1-P)*P_zone/2),
                Exponential(21.0u"yr", (1-P_zone)/2),
                Exponential(10.0u"hr"),
                Exponential(10.0u"hr"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d"),
                Exponential(14.0u"d")
            ]
        )
    else
        error("Invalid process type: $process_type. Use 'SemiMarkov' or 'Markov'.")
    end
    
    return std
end

"""
    generate_component_stds(device_cable_mapping, fault_probabilities, λᶜ; converter_probabilities=nothing, process_type="SemiMarkov", P_zone=1.0)

Generate STDs for all components based on device and cable configurations.

# Arguments
- `device_cable_mapping`: Mapping of devices to cable configurations
- `fault_probabilities`: Fault clearing probabilities for each device/cable combination
- `λᶜ`: Cable failure rate
- `converter_probabilities`: Optional converter fault clearing probabilities
- `process_type`: "SemiMarkov" or "Markov"
- `P_zone`: Zone protection factor

# Returns
- `Dict`: Nested dictionary of STDs organized by device type and cable ID
"""
function generate_component_stds(device_cable_mapping, fault_probabilities, λᶜ; 
                                converter_probabilities=nothing, process_type="SemiMarkov", P_zone=1.0)
    stds = Dict()
    
    for (component_type, device_configs) in device_cable_mapping
        stds[component_type] = Dict()
        
        for (device_type, cable_configs) in device_configs
            stds[component_type][device_type] = Dict()
            
            for (cable_id, cable_lengths) in cable_configs
                P = fault_probabilities[component_type][device_type][cable_id]
                failure_rate = λᶜ * maximum(cable_lengths)
                
                if converter_probabilities !== nothing && haskey(converter_probabilities[component_type][device_type], cable_id)
                    Pc = converter_probabilities[component_type][device_type][cable_id]
                    stds[component_type][device_type][cable_id] = create_component_std(
                        P, Pc, failure_rate; process_type=process_type, P_zone=P_zone
                    )
                else
                    stds[component_type][device_type][cable_id] = create_component_std(
                        P, failure_rate; process_type=process_type, P_zone=P_zone
                    )
                end
            end
        end
    end
    
    return stds
end

"""
    solve_component_stds!(stds, process_class, simulation_time, time_step)

Solve all STDs in the nested dictionary structure.

# Arguments
- `stds`: Nested dictionary of STDs
- `process_class`: Process class (SemiMarkovProcess() or MarkovProcess())
- `simulation_time`: Total simulation time
- `time_step`: Time step for simulation

# Modifies
- `stds`: STDs are solved in-place
"""
function solve_component_stds!(stds, process_class, simulation_time, time_step)
    for (component_type, device_stds) in stds
        for (device_type, cable_stds) in device_stds
            if isa(cable_stds, AbstractDict)
                for (cable_id, std) in cable_stds
                    _MSS.solve!(std, process_class; tsim=simulation_time, dt=time_step)
                end
            else
                _MSS.solve!(cable_stds, process_class; tsim=simulation_time, dt=time_step)
            end
        end
    end
end

"""
    create_reduced_stds(original_stds)

Create reduced state representations from solved STDs.

# Arguments
- `original_stds`: Dictionary of solved STDs

# Returns
- `Dict`: Dictionary of reduced STDs with consolidated states
"""
function create_reduced_stds(original_stds)
    reduced_stds = Dict()
    
    for (component_type, device_stds) in original_stds
        reduced_stds[component_type] = Dict()
        
        for (device_type, cable_stds) in device_stds
            if isa(cable_stds, AbstractDict)
                reduced_stds[component_type][device_type] = Dict()
                for (cable_id, std) in cable_stds
                    reduced_stds[component_type][device_type][cable_id] = create_reduced_std_representation(std)
                end
            else
                reduced_stds[component_type][device_type] = create_reduced_std_representation(cable_stds)
            end
        end
    end
    
    return reduced_stds
end

"""
    create_reduced_std_representation(std_solution)

Create a reduced state representation from a solved STD.

# Arguments
- `std_solution`: Solved STD object

# Returns
- `solvedSTD`: Reduced STD with consolidated states
"""
function create_reduced_std_representation(std_solution)
    time = std_solution.props[:time]
    num_states = 6
    
    # Initialize state probabilities and hazard rates
    state_probs = [zeros(length(time)) for _ in 1:num_states]
    state_hazards = [zeros(length(time)) * unit(std_solution.sprops[1][:h][1]) for _ in 1:num_states]
    
    # Consolidate states based on naming convention
    for n in 1:length(std_solution.sprops)
        state_name = std_solution.sprops[n][:name]
        prob = std_solution.sprops[n][:prob]
        hazard = std_solution.sprops[n][:h]
        
        if occursin("A", state_name)
            # Available state
            state_probs[1] .+= prob
            state_hazards[1] .+= hazard
        elseif occursin("V1", state_name)
            # Voltage fault state
            state_probs[2] .+= vcat(0, prob[2:end])
            state_hazards[2] .+= hazard
        elseif occursin("U2", state_name)
            # Short unavailable state with corrective action
            corrective_prob = _MSS.state_conv(LogNormal(log(2.0)u"hr", 0.25u"hr"), hazard, time, 10000)
            state_probs[3] .+= corrective_prob
            state_hazards[3] .+= hazard
            state_probs[4] .+= vcat(0, prob[2:end]) .- corrective_prob
            state_hazards[4] .+= hazard
        elseif occursin("U3", state_name)
            # Long unavailable state with corrective action
            corrective_prob = _MSS.state_conv(LogNormal(log(20.0)u"hr", 0.05u"hr"), hazard, time, 10000)
            state_probs[5] .+= corrective_prob
            state_hazards[5] .+= hazard
            state_probs[6] .+= vcat(0, prob[2:end]) .- corrective_prob
            state_hazards[6] .+= vcat(0 * unit(std_solution.sprops[1][:h][1]), hazard[1:end-1])
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
