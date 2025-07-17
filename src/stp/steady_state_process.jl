################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models, often found in           #
# reliability engineering.                                                     #
# See https://github.com/timmyfaraday/MultiStateSystems.jl                     #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# types ########################################################################
## abstract types
"""
    SteadyStateProcess <: AbstractMarkovProcess

A steady-state stochastic process type for solving state-transition diagrams
in their steady-state regime. This assumes the system has reached equilibrium
and state probabilities are constant over time.

Use this process type when you're interested in the long-term behavior of
the system rather than its time-dependent evolution.

# Example
```julia-repl
julia> std = STD()
julia> add_states!(std, power = [0u"MW", 100u"MW"], init = [0.0, 1.0])
julia> add_transitions!(std, rate = [0.01u"1/hr", 0.001u"1/hr"], 
                            states = [(2,1), (1,2)])
julia> solve!(std, SteadyStateProcess())
```
"""
mutable struct SteadyStateProcess <: AbstractMarkovProcess end 

# constants ####################################################################
## properties
const steady_state_process_props = [:renewal, :markovian, :steady_state]

# functions ####################################################################
## parameters 
""
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
""
function set_parameters!(std::AbstractSTD, cls::SteadyStateProcess)
    set_rates!(std, cls)
    set_markov_chain_matrix!(std, cls)
end

## stochastic process
"""
    solve!(std::AbstractSTD, cls::SteadyStateProcess; tsim::Number=1.0u"yr", dt::Number=1.0u"d", tol::Real=1e-8)

Solve a state-transition diagram `std` using the steady-state process `cls`.
This computes the steady-state probabilities by solving the system of linear
equations derived from the Markov chain transition matrix.

# Arguments
- `std::AbstractSTD`: The state-transition diagram to solve
- `cls::SteadyStateProcess`: The steady-state process instance
- `tsim::Number`: Simulation time (not used in steady-state, kept for compatibility)
- `dt::Number`: Time step (not used in steady-state, kept for compatibility)
- `tol::Real`: Tolerance for numerical computations, defaults to 1e-8

# Example
```julia-repl
julia> std = STD()
julia> add_states!(std, power = [0u"MW", 100u"MW"], init = [0.0, 1.0])
julia> add_transitions!(std, rate = [0.01u"1/hr", 0.001u"1/hr"], 
                            states = [(2,1), (1,2)])
julia> solve!(std, SteadyStateProcess())
```
"""
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