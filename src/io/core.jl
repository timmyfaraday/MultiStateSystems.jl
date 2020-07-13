################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
#
"""
dim(value::Any) = length(size(value))
function reduce(Prb::Vector,Val::Vector)
    idx = sortperm(Val)
    sPrb, sVal = Prb[idx], Val[idx]
    val = unique(sVal)
    prb = [sum(sPrb[sorted_range(sVal,nv)]) for nv in val]

    return prb, val
end
sorted_range(sVal,nv) = searchsortedfirst(sVal,nv):searchsortedlast(sVal,nv)

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
reduce(kwargs::Iterators.Pairs, nc; exclude=[]) =
    Dict(key => isa(value,Single) ? value : value[nc] for (key,value) in kwargs
                                                      if !in(key,exclude))
