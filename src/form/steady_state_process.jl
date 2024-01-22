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
function set_markov_chain_matrix!(std::AbstractSTD, cls::SteadyStateProcess)
    P = zeros(_MSM.Measurement, ns(std), ns(std))
    for nt in transitions(std)
        rate = get_prop(std, nt, :rate)
        P[Graphs.src(nt),Graphs.dst(nt)] = 
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
function solve!(std::AbstractSTD, cls::SteadyStateProcess;
                tsim::Number=1.0u"yr", dt::Number=1.0u"d", tol::Real=1e-8)
    # set the input
    set_parameters!(std, cls)
    
    # get the input
    t   = [Inf]
    P   = get_prop(std, :P)
    Ns  = ns(std)
    
    # solve the problem
    sol = hcat(P .- _LA.I(Ns),ones(Ns))' \ [zeros(Ns)..., 1.0]
    
    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)
    set_prop!(std, states(std), :prob, sol)
    
    # set the solved status
    set_info!(std, :solved, true)
end