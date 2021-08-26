################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# network (abbr: ntw)
## structs
struct Network{I<:Int} <: AbstractNetwork{I}
    graph::_MG.DiMultigraph{I}

    props::PropDict

    cmp::Vector{PropDict}
    src::Vector{PropDict}
    usr::Vector{PropDict}

    clib::LibDict
    slib::LibDict
    ulib::LibDict
end

## constructors
################################################################################
# WARNING:  The empty constructor needs to be last in order to overwrite the   #
#           empty constructor created by other contructors, see: discourse -   #
#           keyword argument contructor breaks incomplete constructor.         #                                               #
################################################################################
"""
    Network()

An network constructor.

# Example
```julia-repl
julia> ntw = Network()
```
"""
function Network()
    graph = _MG.DiMultigraph(0)
    props = PropDict(:info => NetworkInfo(),
                     :msr  => Set{Symbol}())
    cmp, src, usr = PropDict[], PropDict[], PropDict[]
    clib, slib, ulib = LibDict(), LibDict(), LibDict()

    return Network(graph,props,cmp,src,usr,clib,slib,ulib)
end

## functions
is_directed(::Type{Network}) = true
is_directed(::Type{Network{Int}}) = true
is_directed(ntw::Network) = true

props(ntw::AbstractNetwork) = ntw.props

get_prop(ntw::AbstractNetwork, prop::Symbol) =
    haskey(props(ntw), prop) ? props(ntw)[prop] : ~ ;
has_prop(ntw::AbstractNetwork, prop::Symbol) = haskey(props(ntw), prop)
set_prop!(ntw::AbstractNetwork, prop::Symbol, value) =
    set_props!(ntw, PropDict(prop => value))
set_props!(ntw::AbstractNetwork, dict::PropDict) = merge!(ntw.props, dict)

nv(ntw::AbstractNetwork) = _MG.nv(ntw.graph)
add_vertex!(ntw::AbstractNetwork) = _MG.add_vertex!(ntw.graph)
add_vertex!(ntw::AbstractNetwork, node::Int) =
    !has_vertex(ntw,node) ? add_vertices!(ntw,node-nv(ntw)) : ~ ;
add_vertices!(ntw::AbstractNetwork, n::Int) = _MG.add_vertices!(ntw.graph, n)
has_vertex(ntw::AbstractNetwork, x...) = _MG.has_vertex(ntw.graph, x...)

ne(ntw::AbstractNetwork) = _MG.ne(ntw.graph)
add_edge!(ntw::AbstractNetwork, x...) = _MG.add_edge!(ntw.graph, x...)
has_edge(ntw::AbstractNetwork, x...) = _MG.has_edge(ntw.graph, x...)
mul_edge(ntw::AbstractNetwork,edge::Tuple{Int,Int}) = 
    _MG.mul(ntw.graph, edge[1], edge[2])
mul_edge(graph::_MG.DiMultigraph,edge::Tuple{Int,Int}) = 
    _MG.mul(graph, edge[1], edge[2])
update_lib!(type::Symbol,array::Array,lib::Dict) =
    haskey(lib,array[end][type]) ? push!(lib[array[end][type]],length(array)) :
                                   lib[array[end][type]] = [length(array)] ;

weights(ntw::AbstractNetwork) = _LG.weights(ntw.graph)
has_path(ntw::AbstractNetwork, s_node::Int, u_node::Int) =
    _LG.has_path(ntw.graph, s_node, u_node)
max_paths(graph::_MG.DiMultigraph) = 
    _MG.nv(graph) + _MG.ne(graph, count_mul = true)
nodal_paths(graph::_MG.DiMultigraph, s_node::Int, u_node::Int) =
    sort(_LG.yen_k_shortest_paths(graph, s_node, u_node, 
                             _LG.weights(graph), max_paths(graph)).paths)
mul_path(graph::_MG.DiMultigraph, npath::Array{Int,1}) =
    [1:mul_edge(graph,(npath[ni],npath[ni+1])) for ni in 1:length(npath)-1]

function get_extended_graph(ntw::AbstractNetwork, u_node::Int)
    graph = _MG.copy(ntw.graph)
    _MG.add_vertex!(graph)
    
    for src in ntw.src if has_path(ntw, src[:node], u_node)
        _MG.add_edge!(graph, _MG.nv(graph), src[:node])
    end end
    
    return graph, _MG.nv(graph)
end
function paths(ntw::AbstractNetwork, u_node::Int)
    npaths, cpaths = Vector{Int}[], Vector{Expr}[]

    graph, x_node = get_extended_graph(ntw, u_node)
    for npath in unique!(nodal_paths(graph, x_node, u_node))
        for mul in Iterators.product(mul_path(graph,npath)...)
            cpath = Vector{Expr}()
            # sources
            add_src_expr!(cpath, ntw, npath[2], mul[1])
            # components
            for nn in 2:length(npath)-1
                fr, to, ml = npath[nn], npath[nn+1], mul[nn]
                # fr-node cmp
                add_cmp_expr!(cpath, ntw, fr)
                # edge cmp
                edge = _MG.MultipleEdge(fr, to, ml)
                add_cmp_expr!(cpath, ntw, edge)
                # edge = _MG.MultipleEdge(to, fr, ml) # not necessary in a directed graph
                # add_cmp_expr!(cpath, ntw, edge)
            end
            # end-node cmp
            add_cmp_expr!(cpath, ntw, npath[end])
            
            push!(npaths,npath)
            push!(cpaths,cpath)
    end end
    return npaths, cpaths
end
function ntws(ntw::AbstractNetwork)
    cntr = 1
    ntws = Vector{AbstractNetwork}([ntw])
    while length(ntws) >= cntr
        push!(ntws,setdiff(get_ntw(ntws[cntr]),ntws)...)
        cntr += 1
    end
    for nn in ntws set_msr!(nn) end
    if any(x -> get_info(x,:dependent_sources),ntws)
        ni = findlast(x -> get_info(x,:dependent_sources),ntws)
        ntws[1].props[:source_ugf] = UGF(ntws[ni].props[:msr],ntws[ni].src[1][:std])
        for nn in ntws set_info!(nn,:dependent_sources,true) end
    end
    return reverse(ntws)
end

# info
## structs
mutable struct NetworkInfo{B<:Bool} <: AbstractInfo
    solved::B
    dependent_sources::B
    eval_dep::B

    # default constructor
    NetworkInfo() = new{Bool}(false,false,false)
end
mutable struct ComponentInfo{B<:Bool} <: AbstractInfo
    eval_dep::B
    eval_dep_id::B

    # default constructor 
    ComponentInfo() = new{Bool}(false,false)
end
mutable struct UserInfo{B<:Bool} <: AbstractInfo 
    eval_dep::B
    eval_dep_id::B

    # default constructor
    UserInfo() = new{Bool}(false,false)
end
mutable struct SourceInfo{B<:Bool} <: AbstractInfo
    eval_dep::B
    eval_dep_id::B

    # default constructor
    SourceInfo() = new{Bool}(false,false)
end

## functions
get_info(ntw::AbstractNetwork, info::Symbol) = getproperty(ntw.props[:info],info)
set_info!(ntw::AbstractNetwork, info::Symbol, value::Bool) =
    setproperty!(ntw.props[:info], info, value)
get_info(prt::PropDict, info::Symbol) = getproperty(prt[:info],info)
set_info!(prt::PropDict, info:: Symbol, value::Bool) =
    setproperty!(prt[:info], info, value)
set_info!(info::AbstractInfo, field::Symbol, value::Bool) =
    setproperty!(info, field, value)

# dependence
## source dependence
### functions
init_source_dep(ntw::AbstractNetwork, kwargs::Iterators.Pairs) =
    if haskey(kwargs, :dep)
        set_info!(ntw, :dependent_sources, kwargs[:dep])
    end

## evaluation dependence
### functions
init_eval_dep_src(ntw::AbstractNetwork, kwargs::Iterators.Pairs, Ns::Int) =
    if haskey(kwargs, :eval_dep)
        set_info!(ntw, :eval_dep, kwargs[:eval_dep])
        return length(ntw.src) .+ 1:Ns
    else
        return nothing
    end
init_eval_dep_cmp(ntw::AbstractNetwork, kwargs::Iterators.Pairs, Nc::Int) =
    if haskey(kwargs, :eval_dep)
        set_info!(ntw, :eval_dep, kwargs[:eval_dep])
        return length(ntw.cmp) .+ 1:Nc
    else
        return nothing
    end
function init_eval_dep_cmp(ntw::AbstractNetwork, Nc::Int)
    set_info!(ntw, :eval_dep, true)
    return length(ntw.cmp) .+ 1:Nc
end
init_eval_dep_usr(ntw::AbstractNetwork, kwargs::Iterators.Pairs, Nu::Int) =
    if haskey(kwargs, :eval_dep)
        set_info!(ntw, :eval_dep, kwargs[:eval_dep])
        return length(ntw.usr) .+ 1:Nu
    else
        return nothing
    end
set_eval_dep!(prop_dict::PropDict, ni::Int, eval_dep_ids) =
    if haskey(prop_dict, :eval_dep)
        prop_dict[:eval_dep_ids] = eval_dep_ids
        prop_dict[:eval_dep_id] = ni == 1 ? false : true ;
    end
set_eval_info!(info::AbstractInfo, prop_dict::PropDict) =
    if haskey(prop_dict, :eval_dep)
        for i_key in [:eval_dep,:eval_dep_id]
            set_info!(info, i_key, prop_dict[i_key])
            delete!(prop_dict, i_key)
        end
    end

# elements (abbr: elm)
## functions
elements(ntw::AbstractNetwork) = Iterators.flatten((ntw.cmp,ntw.src))
elm_range(elm::PropDict) = 
    !get_info(elm,:eval_dep_id) ? length(elm[:ugf].val) : 1 ;
get_idx(ntw::AbstractNetwork) = [1:elm_range(elm) for elm in elements(ntw)]
get_idx_itr(ntw::AbstractNetwork) = Iterators.product(get_idx(ntw)...)
get_val(ntw::AbstractNetwork) = [elm[:ugf].val for elm in elements(ntw)]
get_prb(ntw::AbstractNetwork) = [elm[:ugf].prb for elm in elements(ntw)
                                               if  !get_info(elm,:eval_dep_id)]
get_ntw(ntw::AbstractNetwork) =
    [elm[:ntw][1] for elm in elements(ntw) if haskey(elm,:ntw)]

## component (abbr: cmp)
### functions
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

"""
    add_components!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds a single component to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_components!(ntwᵖʷʳ, edge = (1,2),
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
        
        update_lib!(:edge,ntw.cmp,ntw.clib)
    end

    return true
end
function add_component!(ntw::AbstractNetwork, node::Int, dict::Dict=PropDict())
    add_vertex!(ntw, node)

    info = ComponentInfo()
    set_eval_info!(info, dict)

    push!(ntw.cmp, PropDict(:node => node, :info => info, dict...))
    
    update_lib!(:node, ntw.cmp, ntw.clib)
end
function add_component!(ntw::AbstractNetwork, edge::Tuple{Int,Int}, dict::Dict=PropDict())
    add_vertex!(ntw, maximum(edge))
    add_edge!(ntw, edge[1], edge[2])
    edge = _MG.MultipleEdge(edge[1], edge[2], mul_edge(ntw, edge))

    info = ComponentInfo()
    set_eval_info!(info, dict)

    push!(ntw.cmp,Dict(:edge => edge, :info => info, dict...))
    
    update_lib!(:edge, ntw.cmp, ntw.clib)
end
function add_bidirectional_component!(ntw::AbstractNetwork; kwargs...)
    haskey(kwargs, :edge) || return false

    eval_dep_ids = init_eval_dep_cmp(ntw, 2)

    for (ne,edge) in enumerate([kwargs[:edge], reverse(kwargs[:edge])])
        prop_dict = reduce(kwargs, 1, excl=[:edge])
        prop_dict[:eval_dep] = true
        set_eval_dep!(prop_dict, ne, eval_dep_ids)
        add_component!(ntw, edge, prop_dict)
    end

    return true
end
function add_bidirectional_component!(ntw::AbstractNetwork, edge::Tuple{Int,Int}, dict::Dict=PropDict())
    if !haskey(dict,:eval_dep_ids) 
        eval_dep_ids = init_eval_dep_cmp(ntw, 2)
    end
    for (ne,edge) in enumerate([edge, reverse(edge)])
        if !haskey(dict,:eval_dep_ids)
            dict[:eval_dep] = true
            set_eval_dep!(dict, ne, eval_dep_ids)
        end
        add_component!(ntw, edge, dict)
    end
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

            add_component!(ntw, node[ni], prop_dict)
    end end
    if haskey(kwargs,:edge)
        edge, Ne = kwargs[:edge], length(kwargs[:edge])

        eval_dep_ids = init_eval_dep_cmp(ntw, kwargs, Ne)

        for ni in 1:Ne
            prop_dict = reduce(kwargs, ni, excl=[:edge])
            set_eval_dep!(prop_dict, ni, eval_dep_ids)

            add_component!(ntw, edge[ni], prop_dict) 
    end end
    return true
end
function add_bidirectional_components!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs, :edge)) || return false

    edge, Ne = kwargs[:edge], length(kwargs[:edge])

    eval_dep_ids = init_eval_dep_cmp(ntw, kwargs, Ne)

    for ne in 1:Ne
        prop_dict = reduce(kwargs, ne, excl=[:edge])
        set_eval_dep!(prop_dict, ne, eval_dep_ids)

        add_bidirectional_component!(ntw, edge[ne], prop_dict) 
    end
    return true
end

## source (abbr: src)
### functions
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

    info = UserInfo()
    init_source_dep(ntw, kwargs)
    
    push!(ntw.src, Dict(:info => info, kwargs...))
    
    update_lib!(:node, ntw.src, ntw.slib)
    
    return true
end
function add_source!(ntw::AbstractNetwork, node::Int, dict::Dict=PropDict())
    add_vertex!(ntw, node)
    
    info = UserInfo()
    set_eval_info!(info, dict)
    
    push!(ntw.src,Dict(:node => node, :info => info, dict...))
    
    update_lib!(:node, ntw.src, ntw.slib)
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

        add_source!(ntw, node[ni], prop_dict)
    end
    
    return true
end

# user (abbr: usr)
## functions
nu(ntw::AbstractNetwork) = length(ntw.usr)
nu(ntw::AbstractNetwork,u_node::Int) = length(ntw.ulib[u_node])
usr(ntw::AbstractNetwork) = ntw.usr
usr(ntw::AbstractNetwork,u_node::Int) = ntw.usr[ntw.ulib[u_node]]
usr_ids(ntw::AbstractNetwork) = 1:nu(ntw)
usr_ids(ntw::AbstractNetwork,u_node::Int) = ntw.ulib[u_node]
usr_nodes(ntw::AbstractNetwork) = keys(ntw.ulib)
"""
    add_user!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds a single user to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_source!(ntwᵖʷʳ, node = 1,
                           ind  = [:EENS])
```
"""
function add_user!(ntw::AbstractNetwork; kwargs...)
    haskey(kwargs,:node) || return false
    
    add_vertex!(ntw,kwargs[:node])

    info = UserInfo()

    push!(ntw.usr,PropDict(:info => info, kwargs...))
    update_lib!(:node,ntw.usr,ntw.ulib)
    return true
end
function add_user!(ntw::AbstractNetwork, node::Int, dict::Dict=PropDict())
    add_vertex!(ntw,node)

    info = UserInfo()
    set_eval_info!(info, dict)

    push!(ntw.usr,PropDict(:node => node, :info => info, dict...))
    update_lib!(:node,ntw.usr,ntw.ulib)
end
"""
    add_users!(ntw::MultiStateSystems.AbstractNetwork; kwargs...)

Adds multiple users to the network `ntw` and fills their corresponding
`PropDict` with the named arguments `kwargs`. Either an uniform arguments is
given which holds for all components or an array is given whith specific
argument for each component.

# Example
```julia-repl
julia> ntwᵖʷʳ = ntw()
julia> add_sources!(ntwᵖʷʳ, node = [1,5,8],
                            ind  = [:EENS])
```
"""
function add_users!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs,:node)) || return false

    node, Nn = kwargs[:node], length(kwargs[:node])
    
    eval_dep_ids = init_eval_dep_usr(ntw, kwargs, Nn)

    for ni in 1:Nn
        prop_dict = reduce(kwargs, ni, excl=[:node])
        set_eval_dep!(prop_dict, ni, eval_dep_ids)
        
        add_user!(ntw, node[ni], prop_dict)
    end
    return true
end
