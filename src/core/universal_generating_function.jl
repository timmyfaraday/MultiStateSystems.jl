################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Universal Generating Function (UGF)
"""
# structs
abstract type AbstractUGF end
struct UGF <: AbstractUGF
    msr::Symbol
    prb::Vector
    val::Vector
end

# constructors
function UGF(msr::Symbol)
    prb = [1.0]
    val = [(Inf)u"kg/hr"]

    return UGF(msr,prb,val)
end
function UGF(msr::Symbol,std::AbstractSTD)
    prb, val = reduce(get_sprop(std,:prob),get_sprop(std,msr))

    return UGF(msr,prb,val)
end
################################################################################
# WARNING:  The empty constructor needs to be last in order to overwrite the   #
#           empty constructor created by other contructors, see: discourse -   #
#           keyword argument contructor breaks incomplete constructor.         #                                               #
################################################################################
function UGF()
    msr = :msr
    prb = Vector()
    val = Vector()

    return UGF(msr,prb,val)
end

# functions
function cmp_ugf(msr::Symbol,cmp::PropDict)
    if haskey(cmp,:std) return UGF(msr,cmp[:std]) end
    return UGF(msr)
end
function src_ugf(msr::Symbol,src::PropDict)
    if haskey(src,:dep) return UGF(msr,[1.0],[1.0u"MW"]) end
    if haskey(src,:std) return UGF(msr,src[:std]) end
    return UGF(msr)
end
function set_ugf!(ntw::AbstractNetwork)
    msr = ntw.props[:msr]
    for cmp in cmp(ntw) cmp[:ugf] = cmp_ugf(msr,cmp) end
    for src in src(ntw) src[:ugf] = src_ugf(msr,src) end
end

"""
# Probability Function
"""
probability_function(pr::Vector,idx_itr) =
    vec([prod([pr[ne][ni[ne]] for ne in 1:length(pr)])[1] for ni in idx_itr])

"""
# Structure Function
"""
function cmp_structure_function(ntw::AbstractNetwork,s_node::Int,u_node::Int)
    npaths, cpaths = paths(ntw, s_node, u_node)
    while length(npaths[1]) ≠ 0
        has_duplicate_paths(npaths) ? vertical_reduction!(npaths,cpaths) :
                                      horizontal_reduction!(npaths,cpaths) ;
    end
    return cpaths[1][1]
end
function src_structure_function(ntw::AbstractNetwork,s_node::Int)
    nc, src_idx = length(ntw.cmp), src_ids(ntw,s_node)
    return ns(ntw,s_node) == 1 ? src_expr(nc+src_idx[1]) :
                                 par([src_expr(nc+ns) for ns in src_idx]) ;
end
function set_structure_function!(ntw::Network)
    for u_node in usr_nodes(ntw)
        expr = nothing
        for s_node in src_nodes(ntw)
            exrp_src = src_structure_function(ntw,s_node)
            exrp_cmp = cmp_structure_function(ntw,s_node,u_node)
            exrp = :(min($exrp_src,$exrp_cmp))
            expr == nothing ? expr = :($exrp) : expr = :($expr + $exrp) ;
        end
        for nu in ntw.ulib[u_node] ntw.usr[nu][:str] = expr end
    end
end

"""
## Flow
"""
par(args) = :(+($(args...)))
ser(args) = :(min($(args...)))

unique_elements(id,paths::Array) =
    setdiff(paths[id],(paths[ni] for ni in CartesianIndices(paths) if ni≠id)...)
has_unique_elements(id,paths::Array) = !isempty(unique_elements(id,paths))
has_duplicate_paths(paths::Array) = !(length(paths)==length(unique(paths)))
vertical_combination(paths::Array) =
    [par_seg(unique([paths[np][ni] for np in CartesianIndices(paths)]))
                                   for ni in CartesianIndices(paths[1])]
par_seg(segment::Array) = length(segment) > 1 ? par(segment) : segment[1] ;

function vertical_reduction!(npaths::Array,cpaths::Array)
    for unique_npath in unique(npaths)
        idx = broadcast(==,Ref(unique_npath),npaths)
        if sum(idx) > 1
            # clean-up of the component path
            new_cpath = vertical_combination(cpaths[idx])
            deleteat!(cpaths,idx)
            push!(cpaths,new_cpath)
            # clean-up of the nodal path
            deleteat!(npaths,idx)
            push!(npaths,unique_npath)
    end end
end
function horizontal_reduction!(npaths::Array,cpaths::Array)
    for ni in CartesianIndices(cpaths)
        if has_unique_elements(ni,npaths)
            # clean-up of the nodal path
            idx = findall(in(unique_elements(ni,npaths)),npaths[ni])
            deleteat!(npaths[ni],idx)
            # clean-up of the component paths TODO enumerate over seperate sequences larger than one
            idx = findall(in(unique_elements(ni,cpaths)),cpaths[ni])
            cpaths[ni][idx[1]] = ser([cpaths[ni][id] for id in idx])
            deleteat!(cpaths[ni],idx[2:end])
    end end
end

"""
# Solve
"""
function solve!(ntw::AbstractNetwork; type::Symbol=:steady)
    set_msr!(ntw)
    set_ugf!(ntw)

    set_structure_function!(ntw)

    pr, vl = get_prb(ntw,type = type), get_val(ntw)
    idx_itr = Iterators.product(get_idx(ntw)...)
    Prb = probability_function(pr,idx_itr)
    for usr in ntw.usr
        Val = Number[]

        expr = usr[:str]
        exrp = quote function structure_function(idx,val) $expr end end
        eval(exrp)
        
        for ni in idx_itr push!(Val,Base.invokelatest(structure_function,ni,vl)) end

        prb, val = reduce(Prb,Val)

        if haskey(ntw.src[1],:dep)
            ugf = UGF(:flow,ntw.src[1][:std])
            Prb, Val = kron(prb,[pr[end] for pr in ugf.prb]), kron(val,_UF.ustrip.(ugf.val))
            prb, val = reduce(Prb,Val)
        end

        usr[:std] = STD(prob = prb,flow = val)
        usr[:ugf] = UGF(get_msr(ntw),prb,val)
    end
end
