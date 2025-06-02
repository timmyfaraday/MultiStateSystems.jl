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
ns(std::AbstractSTD) = Graphs.nv(std.graph)
states(std::AbstractSTD) = Graphs.vertices(std.graph)
add_vertex!(std::AbstractSTD) = Graphs.add_vertex!(std.graph)
has_vertex(std::AbstractSTD, x...) = Graphs.has_vertex(std.graph, x...)

props(std::AbstractSTD, s::Int) = get(std.sprops, s, PropDict())
get_sprop(std::AbstractSTD, prop::Symbol) =
    [get_prop(std, ns, prop) for ns in states(std)]
get_prop(std::AbstractSTD, ns::Int, prop::Symbol) = props(std, ns)[prop]
has_prop(std::AbstractSTD, ns::Int, prop::Symbol) = haskey(props(std,ns), prop)
set_prop!(std::AbstractSTD, ns::Int, prop::Symbol, value::Any) =
    set_props!(std,ns,Dict{Symbol,Any}(prop => value))
set_prop!(std::AbstractSTD, states::Base.OneTo{Int}, prop::Symbol, value::Any) =
    for ns in states set_props!(std,ns,Dict{Symbol,Any}(prop => value[ns])) end
set_props!(std::AbstractSTD, ns::Int, prop_dict::Dict) =
    haskey(std.sprops,ns) ? merge!(std.sprops[ns],prop_dict) :
                            std.sprops[ns] = prop_dict ;

## ccdf
ccdf(std::AbstractSTD, ns::Int, φ::Quantity, t::Quantity) =
        1.0 - sum(cdf(get_prop(std, Graphs.Edge(ns,nx), :distr), φ, t)
                for nx in Graphs.outneighbors(std.graph, ns); init = 0.0)

## state
"""
    add_state!(std::MultiStateSystems.AbstractSTD; kwargs...)

Adds a single state to the state-transition diagram `std` and fills its
corresponding `PropDict` with the named arguments `kwargs`.

# Example
```julia-repl
julia> stdᵍᵉⁿ = STD()
julia> add_state!(stdᵍᵉⁿ, name  = "normal operation state",
                          power = 100u"MW",
                          init  = 1.0)
```
"""
function add_state!(std::AbstractSTD; kwargs...)
    test(kwargs) || return false

    add_vertex!(std) || return false
    
    set_prop!(std, ns(std), :info, StateInfo())
    set_props!(std, ns(std), Dict(kwargs...))
    
    return true
end
"""
    add_states!(std::MultiStateSystems.AbstractSTD; kwargs...)

Adds multiple states to the state-transition diagram `std` and fills their
corresponding `PropDict` with the named arguments `kwargs`. Either an uniform
argument is given which holds for all states or an array is given with the
specific argument for each state.

# Example
```julia-repl
julia> stdᵍᵉⁿ = STD()
julia> add_states!(stdᵍᵉⁿ, name  = ["normal operation state","failed state"],
                           power = [100.0u"MW",0.0u"MW"],
                           init  = [1.0,0.0],
                           markovian = true)
```
"""
function add_states!(std::AbstractSTD; kwargs...)
    test(kwargs) || return false

    has_msr(kwargs) ? set_prop!(std, :msr, get_msr(kwargs)) : ~ ;               
    for ni in indices_of(kwargs) add_state!(std; reduce(kwargs, ni)...) end
    
    return true
end
