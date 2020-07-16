################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Network
"""
# structs
struct Network{I<:Int} <: AbstractNetwork{I}
    graph::_MG.Multigraph{I,I}

    props::PropDict

    cmp::Vector{PropDict}
    src::Vector{PropDict}
    usr::Vector{PropDict}

    clib::LibDict
    slib::LibDict
    ulib::LibDict
end

# constructors
################################################################################
# WARNING:  The empty constructor needs to be last in order to overwrite the   #
#           empty constructor created by other contructors, see: discourse -   #
#           keyword argument contructor breaks incomplete constructor.         #                                               #
################################################################################
function Network()
    graph = _MG.Multigraph(0)
    props = PropDict(:info => NetworkInfo())
    cmp, src, usr = PropDict[], PropDict[], PropDict[]
    clib, slib, ulib = LibDict(), LibDict(), LibDict()

    return Network(graph,props,cmp,src,usr,clib,slib,ulib)
end

# functions
is_directed(::Type{Network}) = false
is_directed(::Type{Network{Int}}) = false
is_directed(ntw::Network) = false

props(ntw::AbstractNetwork) = ntw.props

get_info(ntw::AbstractNetwork, info::Symbol) = getproperty(ntw.props[:info],info)
set_info!(ntw::AbstractNetwork, info::Symbol, value::Bool) =
    setproperty!(ntw.props[:info],info,value)

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
    ntw.graph.adjmx[edge[1],edge[2]]

get_msr(ntw::AbstractNetwork) = ntw.props[:msr]
set_msr!(ntw::AbstractNetwork) = set_prop!(ntw,:msr,:flow)

update_lib!(type::Symbol,array::Array,lib::Dict) =
    haskey(lib,array[end][type]) ? push!(lib[array[end][type]],length(array)) :
                                   lib[array[end][type]] = [length(array)] ;

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

"""
## Info
"""
# structs
mutable struct NetworkInfo{B<:Bool} <: AbstractInfo
    solved::B
    dependent_sources::B

    # constructor
    NetworkInfo() = new{Bool}(false,false)
end

"""
## Elements
"""
# functions
elements(ntw::AbstractNetwork) = Iterators.flatten((ntw.cmp,ntw.src))

get_idx(ntw::AbstractNetwork) = [1:length(elm[:ugf].prb) for elm in elements(ntw)]
get_idx_itr(ntw::AbstractNetwork) = Iterators.product(get_idx(ntw)...)
get_val(ntw::AbstractNetwork) = [elm[:ugf].val for elm in elements(ntw)]
get_steady_prb(elm::PropDict) =
    dim(elm[:ugf].prb) > 0 ? [prb[end] for prb in elm[:ugf].prb] : elm[:ugf].prb ;
get_prb(ntw::AbstractNetwork; type::Symbol) =
    return type==:steady ? [get_steady_prb(elm) for elm in elements(ntw)] :
                           [elm[:ugf].prb for elm in elements(ntw)] ;
get_ntw(ntw::AbstractNetwork) =
    [elm[:ntw][1] for elm in elements(ntw) if haskey(elm,:ntw)]
"""
### Component
"""
# functions
nc(ntw::AbstractNetwork) = length(ntw.cmp)
nc(ntw::AbstractNetwork,c_key::UIE) = length(ntw.clib[c_key])
cmp(ntw::AbstractNetwork) = ntw.cmp
cmp(ntw::AbstractNetwork,c_key::UIE) = ntw.cmp[ntw.clib[c_key]]
cmp_ids(ntw::AbstractNetwork) = 1:nc(ntw)
cmp_ids(ntw::AbstractNetwork,c_key::UIE) = ntw.clib[c_key]
cmp_expr(nc::Int) = :(val[$nc][idx[$nc]])
cmp_keys(ntw::AbstractNetwork) = keys(ntw.clib)

function add_component!(ntw::AbstractNetwork; kwargs...)
    (haskey(kwargs,:node) || haskey(kwargs,:edge)) || return false
    if haskey(kwargs,:node)
        node = kwargs[:node]
        add_vertex!(ntw,node)
        push!(ntw.cmp,PropDict(kwargs...))
        update_lib!(:node,ntw.cmp,ntw.clib)
    end
    if haskey(kwargs,:edge)
        edge = kwargs[:edge]
        add_vertex!(ntw,maximum(edge))
        add_edge!(ntw,edge[1],edge[2])
        edge = _MG.MultipleEdge(edge[1],edge[2],mul_edge(ntw,edge))
        push!(ntw.cmp,Dict(:edge => edge, reduce(kwargs,1,exclude=[:edge])...))
        update_lib!(:edge,ntw.cmp,ntw.clib)
    end
end
function add_component!(ntw::AbstractNetwork, node::Int, dict::Dict=PropDict())
    add_vertex!(ntw,node)
    push!(ntw.cmp,PropDict(:node => node, dict...))
    update_lib!(:node,ntw.cmp,ntw.clib)
end
function add_component!(ntw::AbstractNetwork, edge::Tuple{Int,Int}, dict::Dict=PropDict())
    add_vertex!(ntw,maximum(edge))
    add_edge!(ntw,edge[1],edge[2])
    edge = _MG.MultipleEdge(edge[1],edge[2],mul_edge(ntw,edge))
    push!(ntw.cmp,Dict(:edge => edge, dict...))
    update_lib!(:edge,ntw.cmp,ntw.clib)
end
function add_components!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && (haskey(kwargs,:node) || haskey(kwargs,:edge))) || return false
    if haskey(kwargs,:node)
        node = kwargs[:node]
        for ni in 1:length(node)
            add_component!(ntw, node[ni], reduce(kwargs,ni,exclude=[:node]))
    end end
    if haskey(kwargs,:edge)
        edge = kwargs[:edge]
        for ni in 1:length(edge)
            add_component!(ntw, edge[ni], reduce(kwargs,ni,exclude=[:edge]))
    end end
end

"""
### Source
"""
# functions
ns(ntw::AbstractNetwork) = length(ntw.src)
ns(ntw::AbstractNetwork,s_node::Int) = length(ntw.slib[s_node])
src(ntw::AbstractNetwork) = ntw.src
src(ntw::AbstractNetwork,s_node::Int) = ntw.src[ntw.slib[s_node]]
src_ids(ntw::AbstractNetwork) = 1:ns(ntw)
src_ids(ntw::AbstractNetwork,s_node::Int) = ntw.slib[s_node]
src_expr(ns::Int) = :(val[$ns][idx[$ns]])
src_nodes(ntw::AbstractNetwork) = keys(ntw.slib)

function add_source!(ntw::AbstractNetwork; node::Int, kwargs...)
    add_vertex!(ntw,node)
    if haskey(kwargs,:dep) set_info!(ntw,:dependent_sources,kwargs[:dep]) end
    push!(ntw.src,Dict(:node => node, kwargs...))
    update_lib!(:node,ntw.src,ntw.slib)
end
function add_source!(ntw::AbstractNetwork, node::Int, dict::Dict=PropDict())
    add_vertex!(ntw,node)
    push!(ntw.src,Dict(:node => node, dict...))
    update_lib!(:node,ntw.src,ntw.slib)
end
function add_sources!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs,:node)) || return false
    node = kwargs[:node]
    if haskey(kwargs,:dep) set_info!(ntw,:dependent_sources,kwargs[:dep]) end
    if dim(node) > 0 for ni in 1:length(node)
        add_source!(ntw, node[ni], reduce(kwargs,ni,exclude=[:node]))
    end
    else for ni in indices_of(kwargs)
        add_source!(ntw, node, reduce(kwargs,ni,exclude=[:node]))
    end end
end

"""
## User
"""
# functions
nu(ntw::AbstractNetwork) = length(ntw.usr)
nu(ntw::AbstractNetwork,u_node::Int) = length(ntw.ulib[u_node])
usr(ntw::AbstractNetwork) = ntw.usr
usr(ntw::AbstractNetwork,u_node::Int) = ntw.usr[ntw.ulib[u_node]]
usr_ids(ntw::AbstractNetwork) = 1:nu(ntw)
usr_ids(ntw::AbstractNetwork,u_node::Int) = ntw.ulib[u_node]
usr_nodes(ntw::AbstractNetwork) = keys(ntw.ulib)

function add_user!(ntw::AbstractNetwork; node::Int, kwargs...)
    add_vertex!(ntw,node)
    push!(ntw.usr,PropDict(:node => node, kwargs...))
    update_lib!(:node,ntw.usr,ntw.ulib)
end
function add_user!(ntw::AbstractNetwork, node::Int, dict::Dict=PropDict())
    add_vertex!(ntw,node)
    push!(ntw.usr,PropDict(:node => node, dict...))
    update_lib!(:node,ntw.usr,ntw.ulib)
end
function add_users!(ntw::AbstractNetwork; kwargs...)
    (test(kwargs) && haskey(kwargs,:node)) || return false
    node = kwargs[:node]
    for ni in 1:length(node)
        add_user!(ntw, node[ni], reduce(kwargs,ni,exclude[:node]))
    end
end

"""
## Pathså
"""
# functions
weights(ntw::AbstractNetwork) = _LG.weights(ntw.graph)
max_paths(ntw::AbstractNetwork) = _MG.nv(ntw.graph) + _MG.ne(ntw.graph,true)
nodal_paths(ntw::AbstractNetwork,s_node::Int,u_node::Int) =
    _LG.yen_k_shortest_paths(_LG.Graph(ntw.graph.adjmx),s_node,u_node,
                             weights(ntw),max_paths(ntw)).paths
mul_path(ntw::AbstractNetwork,npath::Array{Int,1}) =
    [1:mul_edge(ntw,(npath[ni],npath[ni+1])) for ni in 1:length(npath)-1]
function paths(ntw::AbstractNetwork, s_node::Int, u_node::Int)
    npaths, cpaths = Array{Int,1}[], Array{Expr,1}[]
    for npath in nodal_paths(ntw,s_node,u_node)
        for mul in Iterators.product(mul_path(ntw,npath)...)
            cpath = Array{Expr,1}()
            for nn in 1:length(npath)-1
                node = npath[nn]
                haskey(ntw.clib,node) ? push!(cpath,cmp_expr(ntw.clib[node][1])) : ~ ;
                edge = _MG.MultipleEdge(npath[nn],npath[nn+1],mul[nn])
                haskey(ntw.clib,edge) ? push!(cpath,cmp_expr(ntw.clib[edge][1])) : ~ ;
                edge = _MG.MultipleEdge(npath[nn+1],npath[nn],mul[nn])
                haskey(ntw.clib,edge) ? push!(cpath,cmp_expr(ntw.clib[edge][1])) : ~ ;
            end
            node = npath[end]
            haskey(ntw.clib,node) ? push!(cpath,value(ntw.clib[node][1])) : ~ ;
            push!(npaths,npath)
            push!(cpaths,cpath)
    end end
    return npaths, cpaths
end
