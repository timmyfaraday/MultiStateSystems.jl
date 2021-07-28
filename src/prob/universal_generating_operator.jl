################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# Universal Generating Operator
"""
    solve!(ntw::MultiStateSystems.AbstractNetwork)

This function determines the universal generating function related to the output
of all users `usr` of a network `ntw`.
"""
function solve!(ntw::AbstractNetwork)
    for nn in ntws(ntw) solve_network!(nn) end
    for usr in ntw.usr if haskey(usr,:ind)
        for ind in usr[:ind]
            if ind == :GRA usr[:GRA] = GRA(usr) end
    end end end
    if get_info(ntw, :dependent_sources)
        s_ugf = ntw.props[:source_ugf]
        s_prb, s_val = [pr[end] for pr in s_ugf.prb], s_ugf.val
        for usr in ntw.usr
            n_ugf = usr[:ugf]
            n_prb, n_val = [pr[end] for pr in n_ugf.prb], n_ugf.val
            usr[:mat] = s_prb * n_prb'
            val, prb = reduce(kron(n_val,ustrip.(s_val)), kron(n_prb,s_prb))

            msr = Expr(:kw, get_msr(ntw)[1], n_val)
            usr[:std] = eval(:(STD(prob = $(n_prb), $msr)))
            usr[:ugf] = UGF(get_msr(ntw)[1],val,prb)
    end end
    for usr in ntw.usr if haskey(usr,:ind)
        for ind in usr[:ind]
            if ind == :EENS usr[:EENS] = EENS(usr) end
    end end end
end

## Probability Function
# probability_function(pr::Vector,idx_itr) =
#     vec([prod([pr[ne][ni[ne]] for ne in 1:length(pr)])[1] for ni in idx_itr])
probability_function(pr::Vector) = 
    vec([prod(np) for np in Iterators.product(pr...)])

## Structure Function
function user_structure_function(ntw::AbstractNetwork, u_node::Int)
    npaths, cpaths = paths(ntw, u_node)
    cpaths != [[]] || return nothing
    while length(npaths[1]) ≠ 0
        has_duplicate_paths(npaths) ? vertical_reduction!(npaths,cpaths) :
                                      horizontal_reduction!(npaths,cpaths) ;
    end
    return cpaths[1][1]
end
function set_structure_function!(ntw::Network)
    for u_node in usr_nodes(ntw)
        expr = user_structure_function(ntw, u_node)
        for nu in ntw.ulib[u_node] ntw.usr[nu][:str] = :($expr) end
    end
end

### Flow
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
            if length(idx) > 1
                cpaths[ni][idx[1]] = ser([cpaths[ni][id] for id in idx])
            end
            deleteat!(cpaths[ni],idx[2:end])
    end end
end

function solve_network!(ntw::AbstractNetwork)
    set_msr!(ntw)
    set_ugf!(ntw)
    set_structure_function!(ntw)

    msr = get_msr(ntw)[1]
    pr, vl, idx_itr = get_prb(ntw), get_val(ntw), get_idx_itr(ntw)
    Prb = probability_function(pr)
    skip = Int[]
    
    for (nu,usr) in enumerate(ntw.usr) if !in(nu,skip)
        if get_info(usr,:eval_dep)
            Val = zeros(Number, length(Prb), length(usr[:eval_dep_ids]))
            
            for (nc,nu) in enumerate(usr[:eval_dep_ids])
                expr = ntw.usr[nu][:str]
                exrp = quote function structure_function(idx,val) $expr end end
                eval(exrp)

                for (ni,id) in enumerate(idx_itr) 
                    Val[ni,nc] = Base.invokelatest(structure_function,id,vl) 
                end
            end

            rVal, rPrb = reduce(Val, Prb)

            for (nc,nu) in enumerate(usr[:eval_dep_ids])
                ntw.usr[nu][:ugf] = UGF(msr, rVal[:,nc], rPrb, rdc=false)
                ntw.usr[nu][:std] = STD(ntw.usr[nu][:ugf])
            end

            push!(skip, usr[:eval_dep_ids]...)
        else
            Val = zeros(Number, length(Prb))

            expr = usr[:str]
            exrp = quote function structure_function(idx,val) $expr end end
            eval(exrp)

            for (ni,id) in enumerate(idx_itr) 
                Val[ni] = Base.invokelatest(structure_function,id,vl) 
            end

            usr[:ugf] = UGF(msr, Val, Prb)                                      # Currently, only one measure is supported
            usr[:std] = STD(usr[:ugf])
    end end end
    set_info!(ntw, :solved, true)
end