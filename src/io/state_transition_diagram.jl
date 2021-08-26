################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

## State-Transition Diagram
# structs
struct STD{I<:Int} <: AbstractSTD{I}
    graph::_LG.DiGraph{I}

    props::PropDict
    sprops::Dict{I,PropDict}
    tprops::Dict{_LG.Edge{I},PropDict}
end

# constructors
function STD(Ns::Int)
    graph = _LG.DiGraph(Ns)

    props = PropDict(:info => STDInfo())
    sprops = Dict{Int,PropDict}(ns => PropDict() for ns in 1:Ns)
    tprops = Dict{_LG.Edge{Int},PropDict}()

    return STD(graph, props, sprops, tprops)
end

"""
    STD(ugf::MultiStateSystems.UGF)

An state-transition diagram constructor based on an universal generating 
function `ugf`.

Sets `solved` of the STDInfo to true.

# Example
```julia-repl
julia> ugfᵍᵉⁿ = UGF(:power, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7])
julia> stdᵍᵉⁿ = STD(ugfᵍᵉⁿ)
```
"""
function STD(ugf::AbstractUGF)
    msr, val, prb = ugf.msr, ugf.val, ugf.prb
    msr = Expr(:kw, msr, val)
    
    return eval(:(solvedSTD(prob = $(prb), $msr)))
end

"""
    solvedSTD(;prob::Vector, kwargs...)

An state-transition diagram constructor with given state probabilities `prob`
and any number of other arguments `kwargs`.

Sets `solved` of the STDInfo to true.

# Example
```julia-repl
julia> stdᵍᵉⁿ = solvedSTD(prob  = [0.1,0.2,0.7],
                          power = [0.0u"MW",0.0u"MW",2.0u"MW"])
```
"""
function solvedSTD(;prob::Vector, kwargs...)
    isapprox(1.0, sum(prob), rtol=1e-6) || return false

    std = STD(length(prob))
    set_info!(std,:solved,true)
    set_prop!(std,:msr,collect(intersect(MsrSet,keys(kwargs)))[1])
    set_prop!(std,states(std),:prob,prob)
    for ns in states(std) set_props!(std, ns, reduce(kwargs,ns)) end

    return std
end

################################################################################
# WARNING:  The empty constructor needs to be last in order to overwrite the   #
#           empty constructor created by other contructors, see: discourse -   #
#           keyword argument contructor breaks incomplete constructor.         #                                               #
################################################################################
"""
    STD()

An empty state-transition diagram constructor.

# Example
```julia-repl
julia> stdᵍᵉⁿ = STD()
```
"""
function STD()
    graph = _LG.DiGraph()

    props = PropDict(:info => STDInfo())
    sprops = Dict{Int,PropDict}()
    tprops = Dict{_LG.Edge{Int},PropDict}()

    return STD(graph, props, sprops, tprops)
end

# functions
is_directed(::Type{STD}) = true
is_directed(::Type{STD{I}}) where I = true
is_directed(std::STD) = true

props(std::AbstractSTD) = std.props
props(std::AbstractSTD, s::Int) = get(std.sprops, s, PropDict())
props(std::AbstractSTD, t::_LG.Edge) = get(std.tprops, t, PropDict())

get_prop(std::AbstractSTD, prop::Symbol) =
    haskey(props(std), prop) ? props(std)[prop] : ~ ;
has_prop(std::AbstractSTD, prop::Symbol) = haskey(std.props, prop)
set_prop!(std::AbstractSTD, prop::Symbol, value) =
    set_props!(std, Dict(prop => value))
set_props!(std::AbstractSTD, dict::Dict) = merge!(std.props, dict)

## Info
# structs
mutable struct STDInfo{B<:Bool} <: AbstractInfo
    solved::B
    renewal::B
    markovian::B
    time_homogeneous::B

    # constructor
    STDInfo() = new{Bool}(false,true,true,true)
end
mutable struct StateInfo{B<:Bool} <: AbstractInfo
    renewal::B
    trapping::B

    # constructor
    StateInfo() = new{Bool}(true,false)
end
mutable struct TransInfo{B<:Bool} <: AbstractInfo
    renewal::B
    markovian::B
    time_homogeneous::B

    # constructor
    TransInfo() = new{Bool}(true,true,true)
end

# functions
get_info(std::AbstractSTD, info::Symbol) = getproperty(std.props[:info],info)
set_info!(std::AbstractSTD, info::Symbol, value::Bool) =
    setproperty!(std.props[:info],info,value)
get_info(std::AbstractSTD, ns::Int, info::Symbol) =
    getproperty(std.sprops[ns][:info],info)
set_info!(std::AbstractSTD, ns::Int, info::Symbol, value::Bool) =
    setproperty!(std.sprops[ns][:info],info,value)
get_info(std::AbstractSTD, nt::_LG.Edge, info::Symbol) =
    getproperty(std.tprops[nt][:info],info)
set_info!(std::AbstractSTD, nt::_LG.Edge, info::Symbol, value::Bool) =
    setproperty!(std.tprops[nt][:info],info,value)

## State
# functions
ns(std::AbstractSTD) = _LG.nv(std.graph)
states(std::AbstractSTD) = _LG.vertices(std.graph)
add_vertex!(std::AbstractSTD) = _LG.add_vertex!(std.graph)
has_vertex(std::AbstractSTD, x...) = _LG.has_vertex(std.graph, x...)

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
function add_state!(std::AbstractSTD, prop_dict::Dict)
    add_vertex!(std) || return false
    set_prop!(std, ns(std), :info, StateInfo())                                 # TODO auto capture info prop
    set_props!(std, ns(std), prop_dict)
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
    if has_msr(kwargs) set_prop!(std, :msr, get_msr(kwargs)) end                        # TODO auto capture info prop
    for ni in indices_of(kwargs) add_state!(std,reduce(kwargs,ni)) end
    return true
end

## Transition
# functions
nt(std::AbstractSTD) = _LG.ne(std.graph)
transitions(std::AbstractSTD) = _LG.edges(std.graph)
has_edge(std::AbstractSTD, x...) = _LG.has_edge(std.graph, x...)
add_edge!(std::AbstractSTD, x...) = _LG.add_edge!(std.graph, x...)

get_tprop(std::AbstractSTD, prop::Symbol) =
    Dict(nt => get_prop(std,nt,prop) for nt in transitions(std))
get_prop(std::AbstractSTD, nt::_LG.Edge, prop::Symbol) = props(std,nt)[prop]
has_prop(std::AbstractSTD, nt::_LG.Edge, prop::Symbol) =
    haskey(props(std,nt), prop)
set_prop!(std::AbstractSTD, nt::_LG.Edge, prop::Symbol, value::Any) =
    set_props!(std,nt,Dict{Symbol,Any}(prop => value))
set_props!(std::AbstractSTD, nt::_LG.Edge, prop_dict::Dict) =
    haskey(std.tprops,nt) ? merge!(std.tprops[nt],prop_dict) :
                            std.tprops[nt] = prop_dict ;

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
    test(kwargs) || return false
    haskey(kwargs,:states) || return false
    edge = _LG.Edge(kwargs[:states])
    add_edge!(std, edge) || return false
    set_prop!(std, edge, :info, TransInfo())                                    # TODO auto capture info prop
    set_props!(std, edge, Dict(kwargs...))
    return true
end
function add_transition!(std::AbstractSTD, crd::Tuple{Int,Int}, prop_dict::Dict)
    edge = _LG.Edge(crd)
    add_edge!(std, edge) || return false
    set_prop!(std, edge, :info, TransInfo())                                    # TODO auto capture info prop
    set_props!(std, edge, prop_dict)
    # update_info!(std, edge)
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
        add_transition!(std,crd,reduce(kwargs, ni, excl=[:states]))             # TODO auto capture info prop
    end
    return true
end
