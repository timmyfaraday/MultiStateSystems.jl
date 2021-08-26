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
            usr[:std] = eval(:(solvedSTD(prob = $(n_prb), $msr)))
            usr[:ugf] = UGF(get_msr(ntw)[1],val,prb)
    end end
    for usr in ntw.usr if haskey(usr,:ind)
        for ind in usr[:ind]
            if ind == :EENS usr[:EENS] = EENS(usr) end
    end end end
end

## Probability Function
probability_function(pr::Vector) = 
    vec([prod(np) for np in Iterators.product(pr...)])

## Structure Function
function user_structure_function(ntw::AbstractNetwork, u_node::Int)
    npaths, cpaths = paths(ntw, u_node)
    cpaths != [[]] || return nothing
    while length(npaths[1]) ≠ 0
        has_duplicate_paths(npaths) ? parallel_reduction!(npaths, cpaths) : ~ ;
        has_unique_sequence(npaths) ? series_reduction!(npaths, cpaths) : ~ ;
        while has_bridge_sequence(npaths) bridge_reduction!(npaths, cpaths) end
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
#### Parallel
par(args) = :(+($(args...)))
has_duplicate_paths(paths::Vector{<:Vector}) = 
    !(length(paths) == length(unique(paths)))
parallel_combination(paths::Vector{<:Vector}) =
    [par_seg(unique([paths[np][ni] for np in 1:length(paths)]))
                                   for ni in 1:length(paths[1])]
par_seg(segment::Vector) = length(segment) > 1 ? par(segment) : segment[1] ;
function parallel_reduction!(npaths::Vector{<:Vector}, cpaths::Vector{<:Vector})
    for unique_npath in unique(npaths)
        idx = broadcast(==,Ref(unique_npath),npaths)
        if sum(idx) > 1
            # clean-up of the component path
            new_cpath = parallel_combination(cpaths[idx])
            deleteat!(cpaths,idx)
            push!(cpaths,new_cpath)
            # clean-up of the nodal path
            deleteat!(npaths,idx)
            push!(npaths,unique_npath)
    end end
end

#### Series 
ser(args) = :(min($(args...)))
has_unique_sequence(npaths::Vector{<:Vector}) = !isempty(unique_sequence(npaths))
check_sequence(seq::Vector, npath::Vector) =
    length(seq) <= length(npath) && !has_subvector(npath, seq) && 
    has_subvector(npath, seq[2:end-1])
unique_elements(id::Int, idx::Vector, paths::Vector{<:Vector}) =
    setdiff(paths[id],(paths[ni] for ni in 1:length(paths) if !in(ni,idx))...)
function unique_sequence(npaths::Vector{<:Vector})
    length(npaths) > 1 || return [1]

    all_seqs, nogo_seqs = [], []
    for npath in npaths
        for ni in 1:length(npath)-3+1
            seq = npath[ni:ni+3-1]
            push!(all_seqs, seq)
            for xpath in npaths
                if check_sequence(seq, xpath) push!(nogo_seqs, seq) end
            end
        end
    end
    seqs = setdiff(all_seqs, nogo_seqs)
    !isempty(seqs) || return []
    return [np for (np,npath) in enumerate(npaths)
               if  has_subvector(npath,combine_seqs(seqs))]
end
function combine_seqs(seqs) 
    while true 
        new_seqs = []
        for seq in seqs for alt_seq in seqs 
            if seq[2:end] == alt_seq[1:end-1]
                push!(new_seqs,vcat(seq,alt_seq[end]))
            end end 
        end
        !isempty(new_seqs) || break
        seqs = new_seqs
    end
    return seqs[1]
end
function series_reduction!(npaths::Vector{<:Vector}, cpaths::Vector{<:Vector})
    path_idx = unique_sequence(npaths)
    for np in path_idx
        # clean-up nodal paths
        idx = findall(in(unique_elements(np, path_idx, npaths)), npaths[np])
        deleteat!(npaths[np], idx)
        # clean-up of the component paths
        idx = findall(in(unique_elements(np, path_idx, cpaths)), cpaths[np])
        if length(idx) > 1
            cpaths[np][idx[1]] = ser([cpaths[np][id] for id in idx])
        end
        deleteat!(cpaths[np], idx[2:end])
    end
end

#### bridge
delta(arga,argb) = :(-($arga,$argb))
function bridge(paths) 
    left_path = ser(setdiff(paths[1],paths[4]))
    right_path = ser(setdiff(paths[4],paths[1]))

    left_bottom = setdiff(paths[1],paths[2])[1]
    right_bottom = setdiff(paths[4],paths[3])[1]

    temp = setdiff(paths[1],paths[4],[left_bottom])
    left_top = length(temp) == 1 ? temp[1] : ser(temp) ;
    temp = setdiff(paths[4],paths[1],[right_bottom])
    right_top = length(temp) == 1 ? temp[1] : ser(temp) ;
    
    bridge_comp = setdiff(paths[2],union(paths[1],paths[4]))[1]

    delta_left = delta(left_top,left_bottom)
    delta_right = delta(right_top,right_bottom)

    bridge_path = ser([:(abs($delta_left)), :(abs($delta_right)), bridge_comp])
    cond = :(1 * (_UF.ustrip($delta_left) * _UF.ustrip($delta_right) < 0))

    idx_pre = findfirst(x -> x == temp[1], paths[4])
    idx_end = findlast(x -> x == right_bottom, paths[4])

    return vcat(paths[4][1:idx_pre-1],
                par([left_path, right_path, :(*($bridge_path,$cond))]),
                par([left_bottom,right_bottom]),
                paths[4][idx_end+1:end])
end
swapped_paths(npaths) = [(ni,nj) for ni in 1:length(npaths) for nj in ni:length(npaths)
                                 if is_swapped_path(npaths[ni],npaths[nj])]
has_bridge_sequence(npaths::Vector{<:Vector}) = !isempty(bridge_sequence(npaths))
function is_swapped_path(xpath::Vector, ypath::Vector)
    length(xpath) == length(ypath) || return false 
    idx = findall(x -> !iszero(x), xpath .- ypath)
    return  sum(xpath .- ypath) == 0 && 
            sum(.!iszero.(xpath .- ypath)) == 2 && 
            diff(idx)[1] in [1,-1] &&
            sort(xpath[idx]) == sort(ypath[idx])
end
function bridge_sequence(npaths::Vector{<:Vector})
    !isempty(swapped_paths(npaths)) || return []
    (ni,nj) = swapped_paths(npaths)[1]
    long_left_path, long_right_path = npaths[ni], npaths[nj]
    id = findfirst(x -> !iszero(x), npaths[ni] .- npaths[nj])
    short_left_path = long_right_path[1:end .!= id]
    short_right_path = long_left_path[1:end .!= id]
    short_left_path[1:end .!= id] == short_right_path[1:end .!= id] || return []
    selected_paths = [short_left_path, long_left_path, long_right_path, short_right_path]
    idx = [findfirst(x -> x == 1,broadcast(==,Ref(path),npaths)) for path in selected_paths]
    return idx, short_left_path[1:end .!= id]
end
function bridge_reduction!(npaths::Vector{<:Vector}, cpaths::Vector{<:Vector})
    idx, npath = bridge_sequence(npaths)
    # clean-up nodal paths
    push!(npaths, npath)
    deleteat!(npaths, sort(idx))
    # clean-up component paths
    push!(cpaths, bridge(cpaths[idx]))
    deleteat!(cpaths, sort(idx))
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

# ## Structure Function
# function user_structure_function(ntw::AbstractNetwork, u_node::Int)
#     npaths, cpaths = paths(ntw, u_node)
#     cpaths != [[]] || return nothing
#     while length(npaths[1]) ≠ 0
#         has_duplicate_paths(npaths) ? vertical_reduction!(npaths,cpaths) :
#                                       horizontal_reduction!(npaths,cpaths) ;
#     end
#     return cpaths[1][1]
# end
# function set_structure_function!(ntw::Network)
#     for u_node in usr_nodes(ntw)
#         expr = user_structure_function(ntw, u_node)
#         for nu in ntw.ulib[u_node] ntw.usr[nu][:str] = :($expr) end
#     end
# end

# ### Flow
# par(args) = :(+($(args...)))
# ser(args) = :(min($(args...)))

# unique_elements(id,paths::Array) =
#     setdiff(paths[id],(paths[ni] for ni in CartesianIndices(paths) if ni≠id)...)
# has_unique_elements(id,paths::Array) = !isempty(unique_elements(id,paths))
# has_duplicate_paths(paths::Array) = !(length(paths)==length(unique(paths)))
# vertical_combination(paths::Array) =
#     [par_seg(unique([paths[np][ni] for np in CartesianIndices(paths)]))
#                                    for ni in CartesianIndices(paths[1])]
# par_seg(segment::Array) = length(segment) > 1 ? par(segment) : segment[1] ;

# function vertical_reduction!(npaths::Vector, cpaths::Vector)
#     for unique_npath in unique(npaths)
#         idx = broadcast(==,Ref(unique_npath),npaths)
#         if sum(idx) > 1
#             # clean-up of the component path
#             new_cpath = vertical_combination(cpaths[idx])
#             deleteat!(cpaths,idx)
#             push!(cpaths,new_cpath)
#             # clean-up of the nodal path
#             deleteat!(npaths,idx)
#             push!(npaths,unique_npath)
#     end end
# end
# function horizontal_reduction!(npaths::Vector, cpaths::Vector)
#     for ni in CartesianIndices(cpaths)
#         if has_unique_elements(ni,npaths)
#             # clean-up of the nodal path
#             idx = findall(in(unique_elements(ni,npaths)),npaths[ni])
#             deleteat!(npaths[ni],idx)
#             # clean-up of the component paths TODO enumerate over seperate sequences larger than one
#             idx = findall(in(unique_elements(ni,cpaths)),cpaths[ni])
#             if length(idx) > 1
#                 cpaths[ni][idx[1]] = ser([cpaths[ni][id] for id in idx])
#             end
#             deleteat!(cpaths[ni],idx[2:end])
#     end end
# end
# function unique_sequence(paths::Vector{<:Vector})
#     length(paths) > 1 || return [1] 
#     nogo_seqs, all_seqs = [], []
#     for path in paths
#         for ni in 1:length(path)-3+1
#             seq = path[ni:ni+3-1]
#             push!(all_seqs,seq)
#             for npath in paths
#                 if length(npath) >= length(seq) && has_subvector(npath, seq[2:end-1]) && !has_subvector(npath, seq)
#                     push!(nogo_seqs,seq)
#                 end
#             end
#         end
#     end
#     seqs = setdiff(all_seqs,nogo_seqs)
#     res = [[np for (np,path) in enumerate(paths) if has_subvector(path,seq)] for seq in seqs]
#     sort!(res, by=length)
#     return first(res)
# end