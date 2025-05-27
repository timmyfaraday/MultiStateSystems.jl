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
nc(ntw::AbstractNetwork) = length(ntw.cmp)
nc(ntw::AbstractNetwork, c_key::UIE) = length(ntw.clib[c_key])
cmp(ntw::AbstractNetwork) = ntw.cmp
cmp(ntw::AbstractNetwork, c_key::UIE) = ntw.cmp[ntw.clib[c_key][1]]
cmp_id(ntw::AbstractNetwork, c_key::UIE) = ntw.clib[c_key][1]
cmp_ix(ntw::AbstractNetwork, c_key::UIE) = 
    get_info(cmp(ntw, c_key), :eval_dep) ?
        first(cmp(ntw, c_key)[:eval_dep_ids]) :
        cmp_id(ntw, c_key) ;
cmp_ids(ntw::AbstractNetwork) = 1:nc(ntw)
cmp_ids(ntw::AbstractNetwork, c_key::UIE) = ntw.clib[c_key]
cmp_keys(ntw::AbstractNetwork) = keys(ntw.clib)
cmp_expr(nc::Int) = :(val[$nc][idx[$nc]])
cmp_expr(nc::Int, ni::Int) = :(val[$nc][idx[$ni]])
add_cmp_expr!(cpath::Vector{Expr}, ntw::AbstractNetwork, c_key::UIE) =
    haskey(ntw.clib, c_key) ? push!(cpath, cmp_expr(cmp_id(ntw, c_key),
                                                    cmp_ix(ntw, c_key))) : ~ ;

## dependence 
init_eval_dep_cmp(ntw::AbstractNetwork, kwargs::Iterators.Pairs, Nc::Int) =
    if haskey(kwargs, :eval_dep)
        set_info!(ntw, :eval_dep, kwargs[:eval_dep])
        return length(ntw.cmp) .+ (1:Nc)
    else
        return nothing
    end
""
function init_eval_dep_cmp(ntw::AbstractNetwork, Nc::Int)
    set_info!(ntw, :eval_dep, true)
    return length(ntw.cmp) .+ (1:Nc)
end

## component (directional)
"""
    add_component!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds a single component to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_component!(ntwᵖʷʳ, edge = (1,2),
                              name = "cable 1",
                              std  = STD(power = [0u"MW",1500u"MW"],
                                         prob  = [0.2,0.8]))
```
"""
function add_component!(ntw::AbstractNetwork; kwargs...)
    (haskey(kwargs, :node) || haskey(kwargs, :edge)) || return false

    if haskey(kwargs,:node)
        add_vertex!(ntw, kwargs[:node])

        info = ComponentInfo()

        push!(ntw.cmp, PropDict(:info => info, kwargs...))

        update_lib!(:node, ntw.cmp, ntw.clib)
    end
    
    if haskey(kwargs,:edge)
        edge = kwargs[:edge]
        add_vertex!(ntw, maximum(edge))
        add_edge!(ntw, edge[1], edge[2])
        edge = _MG.MultipleEdge(edge[1], edge[2], mul_edge(ntw,edge))

        info = ComponentInfo()

        push!(ntw.cmp, Dict(:edge => edge, :info => info, reduce(kwargs, 1, excl=[:edge])...))
        
        update_lib!(:edge, ntw.cmp, ntw.clib)
    end

    return true
end
"""
    add_components!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds multiple components to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`. Either an uniform arguments is
given which holds for all components or an array is given whith specific
argument for each component.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_components!(ntwᵖʷʳ, edge = [(1,2),(1,2),(2,3)],
                               name = ["cable 1","cable 2","cable 3"],
                               std  = [STD(power = [0u"MW",1500u"MW"],
                                           prob  = [0.2,0.8]),
                                       STD(power = [0u"MW",2000u"MW"],
                                           prob  = [0.4,0.6]),
                                       STD(power = [0u"MW",1800u"MW",4000u"MW"],
                                           prob = [0.1,0.2,0.7])])
```
"""
function add_components!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && (haskey(kwargs, :node) || haskey(kwargs, :edge))) || return false

    if haskey(kwargs, :node)
        node, Nn = kwargs[:node], length(kwargs[:node])

        eval_dep_ids = init_eval_dep_cmp(ntw, kwargs, Nn)

        for ni in 1:Nn
            prop_dict = reduce(kwargs, ni, excl=[:node])
            set_eval_dep!(prop_dict, ni, eval_dep_ids)

            add_component!(ntw; node=node[ni], prop_dict...)
    end end

    if haskey(kwargs,:edge)
        edge, Ne = kwargs[:edge], length(kwargs[:edge])

        eval_dep_ids = init_eval_dep_cmp(ntw, kwargs, Ne)

        for ni in 1:Ne
            prop_dict = reduce(kwargs, ni, excl=[:edge])
            set_eval_dep!(prop_dict, ni, eval_dep_ids)

            add_component!(ntw; edge=edge[ni], prop_dict...) 
    end end
    return true
end

## component (bidirectional)
""
function add_bidirectional_component!(ntw::AbstractNetwork; kwargs...)
    haskey(kwargs, :edge) || return false

    eval_dep_ids = init_eval_dep_cmp(ntw, 2)

    for (ne,edge) in enumerate([kwargs[:edge], reverse(kwargs[:edge])])
        prop_dict = reduce(kwargs, 1, excl=[:edge])
        prop_dict[:eval_dep] = true
        set_eval_dep!(prop_dict, ne, eval_dep_ids)

        add_component!(ntw; edge=edge, prop_dict...)
    end

    return true
end
""
function add_bidirectional_components!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs, :edge)) || return false

    edge, Ne = kwargs[:edge], length(kwargs[:edge])

    eval_dep_ids = init_eval_dep_cmp(ntw, kwargs, Ne)

    for ne in 1:Ne
        prop_dict = reduce(kwargs, ne, excl=[:edge])
        set_eval_dep!(prop_dict, ne, eval_dep_ids)

        add_bidirectional_component!(ntw; edge=edge[ne], prop_dict...) 
    end

    return true
end
