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

# functions ####################################################################
## general
nt(std::AbstractSTD) = Graphs.ne(std.graph)
transitions(std::AbstractSTD) = Graphs.edges(std.graph)
has_edge(std::AbstractSTD, x...) = Graphs.has_edge(std.graph, x...)
add_edge!(std::AbstractSTD, x...) = Graphs.add_edge!(std.graph, x...)

props(std::AbstractSTD, t::Graphs.Edge) = get(std.tprops, t, PropDict())
get_tprop(std::AbstractSTD, prop::Symbol) =
    Dict(nt => get_prop(std,nt,prop) for nt in transitions(std))
get_prop(std::AbstractSTD, nt::Graphs.Edge, prop::Symbol) = props(std,nt)[prop]
has_prop(std::AbstractSTD, nt::Graphs.Edge, prop::Symbol) =
    haskey(props(std,nt), prop)
set_prop!(std::AbstractSTD, nt::Graphs.Edge, prop::Symbol, value::Any) =
    set_props!(std,nt,Dict{Symbol,Any}(prop => value))
set_props!(std::AbstractSTD, nt::Graphs.Edge, prop_dict::Dict) =
    haskey(std.tprops,nt) ? merge!(std.tprops[nt],prop_dict) :
                            std.tprops[nt] = prop_dict ;

## transition
"""
    add_transition!(std::MultiStateSystems.AbstractSTD; kwargs...)

Adds a single transitions to the state-transition diagram `std` and fills its
corresponding `PropDict` with the named arguments `kwargs`. One obligatory named
argument is `:states`, describing the tuple (fr,to) of the from- and to-state.

# Example
```julia-repl
julia> stdᵍᵉⁿ = STD()
julia> add_states!(stdᵍᵉⁿ, name  = ["normal operation state","failed state"],
                           power = [100.0u"MW",0.0u"MW"],
                           init  = [1.0,0.0],
                           markovian = true)
julia> add_transition!(stdᵍᵉⁿ, rate = 0.001u"1/hr",
                               states = (1,2))
```
"""
function add_transition!(std::AbstractSTD; kwargs...)
    test(kwargs) && haskey(kwargs,:states) || return false

    edge = Graphs.Edge(kwargs[:states])

    add_edge!(std, edge) || return false
    set_prop!(std, edge, :info, TransInfo())
    set_props!(std, edge, Dict(kwargs...))
    
    return true
end
"""
    add_transitions!(std::MultiStateSystems.AbstractSTD; kwargs...)

Adds multiple transitions to the state-transition diagram `std` and fills their
corresponding `PropDict` with the named arguments `kwargs`. Either an uniform
argument is given which holds for all transitions or an array is given with the
specific argument for each transition.

# Example
```julia-repl
julia> stdᵍᵉⁿ = STD()
julia> add_states!(stdᵍᵉⁿ, name  = ["normal operation state","failed state"],
                           power = [100.0u"MW",0.0u"MW"],
                           init  = [1.0,0.0],
                           markovian = true)
julia> add_transitions!(stdᵍᵉⁿ, rate = [0.001u"1/hr",0.01u"1/hr"],
                                states = [(1,2),(2,1)])
```
!!! note
    If the `:states` argument is not provided in the `add_transitions!`
    function, the from- and to-states will be determined based on the other
    arguments.
# Example (Alternative)
```julia-repl
julia> add_transitions!(stdᵍᵉⁿ, rate = [0.000u"1/hr" 0.010u"1/hr"
                                        0.001u"1/hr" 0.000u"1/hr"])
```
"""
function add_transitions!(std::AbstractSTD; kwargs...)
    test(kwargs) || return false

    for ni in indices_of(kwargs)
        crd = haskey(kwargs,:states) ? kwargs[:states][ni] : (ni[1],ni[2]) ;
        
        add_transition!(std, crd, reduce(kwargs, ni, excl=[:states])...) 
    end
    
    return true
end
