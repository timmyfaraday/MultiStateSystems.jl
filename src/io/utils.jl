################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Arrays
"""
dim(value::Any) = length(size(value))
function reduce(Val::Vector,Prb::Vector)
    idx = sortperm(Val)
    sPrb, sVal = Prb[idx], Val[idx]
    val = unique(sVal)
    prb = [sum(sPrb[sorted_range(sVal,nv)]) for nv in val]

    return val, prb
end
function reduce(Val::Array, Prb::Vector)
    idx = sortperm(collect(eachrow(Val)))
    sVal, sPrb = Val[idx,:], Prb[idx]
    val = unique(sVal, dims=1)
    cVal = collect(eachrow(sVal))
    prb = [sum(sPrb[sorted_range(cVal,nv)]) for nv in eachrow(val)]

    return val, prb
end
sorted_range(sVal,nv) = searchsortedfirst(sVal,nv):searchsortedlast(sVal,nv)

function update!(v1::Vector{Float64},v2::Vector{Float64})
    length(v1) > length(v2) ? push!(v2, zeros(length(v1) - length(v2))...) :
                              push!(v1, zeros(length(v2) - length(v1))...) ;
    v1 .+= v2
end
function shift!(array::Vector{Float64},step::Number,offset::Number)
    push!(array,0.0)
    array *= step / (step + offset) 
    array[2:end] += offset / (step + offset) * array[1:end-1]
end
"""
# KWARGS
"""
# functions
isa_svr(kwargs::Iterators.Pairs) =
    prod([isa(nv,Single) || isa(nv,Vector) || isa(nv,UnitRange)
        for nv in values(kwargs)])
isa_sm(kwargs::Iterators.Pairs) =
    prod([isa(nv,Single) || isa(nv,Matrix) for nv in values(kwargs)])
function equal_size(kwargs::Iterators.Pairs)
    values_length = [size(value) for value in values(kwargs)
                                 if !isa(value,Single)]
    empty = isempty(values_length)
    equal = all(x->x==first(values_length),values_length)
    # test
    empty || equal || return false
end
function test(kwargs::Iterators.Pairs)
    (isa_svr(kwargs) || isa_sm(kwargs)) || return false                         # TODO Memento logger
    equal_size(kwargs) || return false                                          # TODO Memento logger
end

indices_of(kwargs::Iterators.Pairs) =
    CartesianIndices(kwargs[findfirst(x->!isa(x,Single),values(kwargs))])
reduce(kwargs::Iterators.Pairs, nc; excl=[]) =
    PropDict(key => isa(value,Single) ? value : value[nc] for (key,value) in kwargs
                                                          if !in(key,excl))

"""
# UNITS
"""
Base.log(a::Quantity{T,D,U}) where T <: Real where D where U = log(a.val)unit(a)
Base.exp(a::Quantity{T,D,U}) where T <: Real where D where U = exp(a.val)unit(a)