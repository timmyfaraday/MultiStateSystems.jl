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

# constructors #################################################################
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

    return Network( graph, 
                    props, 
                    cmp, src, usr, 
                    clib, slib, ulib)
end

# functions ####################################################################
## ntw
""
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

## graph
### vertices
nv(ntw::AbstractNetwork) = _MG.nv(ntw.graph)
add_vertex!(ntw::AbstractNetwork) = _MG.add_vertex!(ntw.graph)
add_vertex!(ntw::AbstractNetwork, node::Int) =
    !has_vertex(ntw,node) ? add_vertices!(ntw,node-nv(ntw)) : ~ ;
add_vertices!(ntw::AbstractNetwork, n::Int) = _MG.add_vertices!(ntw.graph, n)
has_vertex(ntw::AbstractNetwork, x...) = _MG.has_vertex(ntw.graph, x...)

### edges
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

### paths
has_path(ntw::AbstractNetwork, s_node::Int, u_node::Int) =
    Graphs.has_path(ntw.graph, s_node, u_node)
max_paths(graph::_MG.DiMultigraph) = 
    _MG.nv(graph) + _MG.ne(graph, count_mul = true)
nodal_paths(graph::_MG.DiMultigraph, s_node::Int, u_node::Int) =
    sort(Graphs.yen_k_shortest_paths(graph, s_node, u_node, 
                             Graphs.weights(graph), max_paths(graph)).paths)
mul_path(graph::_MG.DiMultigraph, npath::Array{Int,1}) =
    [1:mul_edge(graph,(npath[ni],npath[ni+1])) for ni in 1:length(npath)-1] 
""
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

### extended graph
function get_extended_graph(ntw::AbstractNetwork, u_node::Int)
    graph = _MG.copy(ntw.graph)
    _MG.add_vertex!(graph)
    
    for src in ntw.src if has_path(ntw, src[:node], u_node)
        _MG.add_edge!(graph, _MG.nv(graph), src[:node])
    end end
    
    return graph, _MG.nv(graph)
end

## elements
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