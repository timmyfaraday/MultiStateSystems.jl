################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# Stochastic Process
"""
    solve!(std::MultiStateSystems.AbstractSTD, tsim::Number; alg::Symbol=:nothing)

This function determines the state probabilities of the state-transition diagram
`std`.

The following optional arguments may be passed, with their respective defaults:
- alg
- dyn
- tsim [1.0u"yr"]
Optionally, the prefered stochastic process may be provided through the named
argument `alg`, otherwise the appropriate stochastic process is determined using
the properties of the state-transition diagram.

# Example
```julia-repl
julia> solve!(stdᵍᵉⁿ, 1000u"hr", alg = :markov_process)
```
"""
function solve!(std::AbstractSTD; alg::Symbol=:markov_process, tsim::Number=1.0u"yr")
    if is_a_steady_state_process(std, alg)  solve_steady_state_process(std)
    if is_a_markov_process(std, alg)        solve_markov_process!(std, tsim) end
    set_info!(std,:solved,true)
end

