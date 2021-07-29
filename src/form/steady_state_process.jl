################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# abstract types
mutable struct SteadyStateProcess <: AbstractMarkovProcess end 

# properties
const steady_state_process_props = [:renewal, :markovian, :steady_state]

# parameters 
function set_rates!(std::AbstractSTD, cls::AbstractMarkovProcess)
    for nt in transitions(std)
        if !has_prop(std, nt, :rate) && has_prop(std, nt, :distr)
            if get_prop(std, nt, :distr) isa Exponential
                set_prop!(std, nt, :rate, rate(get_prop(std, nt, :distr)))
    end end end
end
function set_markov_chain_matrix!(std::AbstractSTD, cls::SteadyStateProcess)
    P = zeros(_MSM.Measurement, ns(std), ns(std))
    for nt in transitions(std)
        rate = get_prop(std, nt, :rate)
        P[_LG.src(nt),_LG.dst(nt)] = 
            ifelse(_UF.unit(rate)==_UF.NoUnits, rate, rate * 1.0u"s" |> u"s/s")
    end
    for ns in 1:ns(std)
        P[ns,ns] = 1.0 - sum(P[ns,:])
    end
    set_prop!(std, :P, ifelse(all(_MSM.uncertainty.(P).== 0.0_UF.unit(first(P))),
                              _MSM.value.(P),
                              P))
end
function set_parameters!(std::AbstractSTD, cls::SteadyStateProcess)
    set_rates!(std, cls)
    set_markov_chain_matrix!(std, cls)
end

# stochastic process
solve_steady_state(std) =
    hcat(get_prop(std, :P) .- _LA.I(ns(std)),ones(ns(std)))' \ [zeros(ns(std))..., 1.0]
function solve!(std::AbstractSTD, cls::SteadyStateProcess)
    set_parameters!(std, cls)

    set_prop!(std, :cls, cls)
    set_prop!(std, :time, [Inf])
    set_prop!(std, states(std), :prob, solve_steady_state(std))

    set_info!(std, :solved, true)
end