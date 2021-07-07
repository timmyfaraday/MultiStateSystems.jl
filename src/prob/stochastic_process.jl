################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# stochastic process
"""
    solve!(std::MultiStateSystems.AbstractSTD; tsim::Number=1.0u"yr", tol::Real=1e-8)

This function determines the state probabilities of the state-transition diagram
`std`. The appropriate stochastic process is determined using the properties of 
the state-transition diagram.


The following optional arguments may be passed, with their respective defaults:
- tsim: simulation duration [1.0u"yr"]
- tol: numerical tolerance [1e-8]

# Example
```julia-repl
julia> solve!(stdᵍᵉⁿ, tsim = 1000u"hr", tol = 1e-6)
```
"""
function solve!(std::AbstractSTD; tsim::Number=1.0u"yr", tol::Real=1e-8)
    if is_a_steady_state_process(std)   
        solve!(std, SteadyStateProcess()) 
    end
    if is_a_markov_process(std) 
        solve!(std, MarkovProcess(), tsim = tsim, tol = tol) 
    end
end

