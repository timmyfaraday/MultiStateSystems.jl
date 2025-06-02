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

# structs ######################################################################
""
struct STD{I<:Int} <: AbstractSTD{I}
    graph::Graphs.DiGraph{I}

    props::PropDict
    sprops::Dict{I,PropDict}
    tprops::Dict{Graphs.Edge{I},PropDict}
end

# constructors #################################################################
""
function STD(Ns::Int)
    graph = Graphs.DiGraph(Ns)

    props = PropDict(:info => STDInfo())
    sprops = Dict{Int,PropDict}(ns => PropDict() for ns in 1:Ns)
    tprops = Dict{Graphs.Edge{Int},PropDict}()

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
function solvedSTD(;prob::Vector, time::Vector=[(Inf)u"yr"], kwargs...)
    length(first(prob)) == length(time) || return false
    all(isapprox.(1.0, sum(prob), rtol=1e-6)) || return false

    std = STD(length(prob))
    set_prop!(std,:time,time)
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
    graph = Graphs.DiGraph()

    props = PropDict(:info => STDInfo())
    sprops = Dict{Int,PropDict}()
    tprops = Dict{Graphs.Edge{Int},PropDict}()

    return STD(graph, props, sprops, tprops)
end

# functions ####################################################################
## general
props(std::AbstractSTD) = std.props

get_prop(std::AbstractSTD, prop::Symbol) =
    haskey(props(std), prop) ? props(std)[prop] : ~ ;
has_prop(std::AbstractSTD, prop::Symbol) = haskey(std.props, prop)
set_prop!(std::AbstractSTD, prop::Symbol, value) =
    set_props!(std, Dict(prop => value))
set_props!(std::AbstractSTD, dict::Dict) = merge!(std.props, dict)