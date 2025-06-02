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
ns(ntw::AbstractNetwork) = length(ntw.src)
ns(ntw::AbstractNetwork,s_node::Int) = length(ntw.slib[s_node])
src(ntw::AbstractNetwork) = ntw.src
src(ntw::AbstractNetwork,s_node::Int, mul::Int) = ntw.src[ntw.slib[s_node][mul]]
src_id(ntw::AbstractNetwork, s_node::Int, mul::Int) = 
    nc(ntw) + ntw.slib[s_node][mul]
src_ix(ntw::AbstractNetwork, s_node::Int, mul::Int) =
    get_info(src(ntw, s_node, mul), :eval_dep) ?
        nc(ntw) + first(src(ntw, s_node, mul)[:eval_dep_ids]) :
        nc(ntw) + ntw.slib[s_node][mul] ;
src_ids(ntw::AbstractNetwork, s_node::Int) = ntw.slib[s_node]
src_ids(ntw::AbstractNetwork) = 1:ns(ntw)
src_expr(ns::Int) = :(val[$ns][idx[$ns]])
src_expr(ns::Int, ni::Int) = :(val[$ns][idx[$ni]])
src_nodes(ntw::AbstractNetwork) = keys(ntw.slib)
add_src_expr!(cpath::Vector{Expr}, ntw::AbstractNetwork, node::Int, mul::Int) =
    haskey(ntw.slib, node) ? push!(cpath, src_expr( src_id(ntw, node, mul),
                                                    src_ix(ntw, node, mul))) : ~ ;

## dependence
init_source_dep(ntw::AbstractNetwork, kwargs::Iterators.Pairs) =
    if haskey(kwargs, :dep)
        set_info!(ntw, :dependent_sources, kwargs[:dep])
    end
init_eval_dep_src(ntw::AbstractNetwork, kwargs::Iterators.Pairs, Ns::Int) =
    if haskey(kwargs, :eval_dep)
        set_info!(ntw, :eval_dep, kwargs[:eval_dep])
        return length(ntw.src) .+ (1:Ns)
    else
        return nothing
    end

## source
"""
    add_source!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds a single source to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> stdᵍᵉⁿ = solvedSTD(prob = [0.1,0.2,0.7],
                          flow = [0.0u"MW",0.5u"MW",2.0u"MW"])
julia> add_source!(ntwᵖʷʳ, node = 1,
                           name = "generator 1",
                           std  = stdᵍᵉⁿ)
```
"""
function add_source!(ntw::AbstractNetwork; kwargs...)
    haskey(kwargs, :node) || return false

    add_vertex!(ntw, kwargs[:node])

    info = SourceInfo()
    init_source_dep(ntw, kwargs)
    set_eval_info!(info; kwargs...)
    
    push!(ntw.src, Dict(:info => info, kwargs...))
    
    update_lib!(:node, ntw.src, ntw.slib)
    
    return true
end
"""
    add_sources!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds multiple sources to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`. Either an uniform arguments is
given which holds for all components or an array is given whith specific
argument for each component.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> stdᵍᵉⁿ = solvedSTD(prob = [0.1,0.2,0.7],
                          flow = [0.0u"MW",0.5u"MW",2.0u"MW"])
julia> add_sources!(ntwᵖʷʳ, node = 1:5,
                            std  = stdᵍᵉⁿ,
                            dep  = true)
```
"""
function add_sources!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs,:node)) || return false

    node, Nn = kwargs[:node], length(kwargs[:node])
    
    init_source_dep(ntw, kwargs)
    eval_dep_ids = init_eval_dep_src(ntw, kwargs, Nn)

    for ni in 1:Nn
        prop_dict = reduce(kwargs, ni, excl=[:node])
        set_eval_dep!(prop_dict, ni, eval_dep_ids)

        add_source!(ntw; node=node[ni], prop_dict...)
    end
    
    return true
end