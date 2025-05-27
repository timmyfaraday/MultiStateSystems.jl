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
## arrays
has_subvector(vec::Vector, sub::Vector) = 
  any([sub == vec[i:(i+length(sub)-1)] for i=1:(length(vec)-length(sub)+1)])
""
function reduce(Val::Vector, Prb::Vector)
    idx = sortperm(Val)
    sPrb, sVal = Prb[idx], Val[idx]
    val = unique(sVal)
    prb = [sum(sPrb[sorted_range(sVal,nv)]) for nv in val]

    return val, prb
end
""
function reduce(Val::Array, Prb::Vector)
    idx = sortperm(collect(eachrow(Val)))
    sVal, sPrb = Val[idx,:], Prb[idx]
    val = unique(sVal, dims=1)
    cVal = collect(eachrow(sVal))
    prb = [sum(sPrb[sorted_range(cVal,nv)]) for nv in eachrow(val)]

    return val, prb
end

## dependence
set_eval_dep!(prop_dict::PropDict, ni::Int, eval_dep_ids) =
    if haskey(prop_dict, :eval_dep)
        prop_dict[:eval_dep_ids] = eval_dep_ids
        prop_dict[:eval_dep_id] = ni == 1 ? false : true ;
    end

## dst
eval(ω::Number,t::Number) = ω
eval(ω::Function,t::Number) = ω(t) |> u"s/s" 

## integral 
weights(N::Int, p::Int) = p == 1 || p == N ? 0.5 : 1.0 ;

## kwargs
indices_of(kwargs::Iterators.Pairs) =
    CartesianIndices(kwargs[findfirst(x->!isa(x,Single),values(kwargs))])
isa_svr(kwargs::Iterators.Pairs) =
    prod([isa(nv,Single) || isa(nv,Vector) || isa(nv,UnitRange)
        for nv in values(kwargs)])
isa_sm(kwargs::Iterators.Pairs) =
    prod([isa(nv,Single) || isa(nv,Matrix) for nv in values(kwargs)])
reduce(kwargs::Iterators.Pairs, nc; excl=[]) =
    PropDict(key => isa(value,Single) ? value : value[nc] for (key,value) in kwargs
                                                          if !in(key,excl))
""
function equal_size(kwargs::Iterators.Pairs)
    values_length = [size(value) for value in values(kwargs)
                                 if !isa(value,Single)]
    empty = isempty(values_length)
    equal = all(x->x==first(values_length),values_length)
    
    empty || equal || return false
end
""
function test(kwargs::Iterators.Pairs)
    (isa_svr(kwargs) || isa_sm(kwargs)) || return false
    equal_size(kwargs) || return false
end

## linear algebra
_LA.:\(A::Matrix{<:Number}, b::Vector{<:Number}) = (ustrip.(A) \ ustrip(b)) * (elunit(b) / elunit(A))

## range
sorted_range(sVal,nv) = searchsortedfirst(sVal,nv):searchsortedlast(sVal,nv)

## units
elunit(A::Matrix{U}) where U = _UF.unit(U)
elunit(b::Vector{U}) where U = _UF.unit(U)
Base.log(a::Quantity{T,D,U}) where T <: Real where D where U = log(a.val)unit(a)
Base.exp(a::Quantity{T,D,U}) where T <: Real where D where U = exp(a.val)unit(a)