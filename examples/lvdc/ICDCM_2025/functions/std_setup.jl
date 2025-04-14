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