################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

## Universal Generating Function
# structs
"""
    UGF

An ugf is a struct containing: a measure `msr`, corresponding values `val` and 
associated probabilities `prb`, with default constructor: 
    
    UGF(msr::Symbol, val::Vector, prb::Vector; red::Bool=true)

A UGF constructor for a specific measure `msr` based on a given value vector 
`val` and associated probability vector `prb`.

This function automatically reduces the state-space to where only unique values
and associated probabilites remain. This behavior is circumvented by passing the
optional argument `rdc=false`.

# Example
```julia-repl
julia> ugfᵍᵉⁿ = UGF(:power, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7])
julia> isequal(ugfᵍᵉⁿ.val, [0.0u"MW",2.0u"MW"])
true
julia> ugfᵍᵉⁿ = UGF(:power, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7], rdc=false)
julia> isequal(ugfᵍᵉⁿ.val, [0.0u"MW",0.0u"MW",2.0u"MW"])
true
```
"""
struct UGF <: AbstractUGF
    msr::Symbol
    val::Vector
    prb::Vector

    # default constructor
    UGF(msr::Symbol, val::Vector, prb::Vector; rdc::Bool=true) =
        rdc ? new(msr,reduce(val,prb)...) : new(msr,val,prb)
end

# constructors
function UGF(msr::Symbol)
    prb = [1.0]
    val = [get_max(msr)]

    return UGF(msr, val, prb)
end

"""
    UGF(msr::Symbol, std::MultiStateSystems.AbstractSTD)

A UGF constructor for a specific measure `msr` based on a given state-transition
diagram `std`.

This function automatically reduces the state-space to where only unique values
and associated probabilites remain.

# Example
```julia-repl
julia> stdᵍᵉⁿ = solvedSTD(prob  = [0.1,0.2,0.7],
                          power = [0.0u"MW",0.0u"MW",2.0u"MW"])
julia> ugfᵍᵉⁿ = UGF(:power, stdᵍᵉⁿ)
julia> isequal(ugfᵍᵉⁿ.val,[0.0u"MW",2.0u"MW"])
true
```
"""
UGF(msr::Symbol, std::AbstractSTD) = 
    UGF(msr,get_sprop(std,msr),get_sprop(std,:prob))
################################################################################
# WARNING:  The empty constructor needs to be last in order to overwrite the   #
#           empty constructor created by other contructors, see: discourse -   #
#           keyword argument contructor breaks incomplete constructor.         #                                               #
################################################################################
function UGF()
    msr = :msr
    prb = Vector()
    val = Vector()

    return UGF(msr,val,prb)
end

# functions
function cmp_ugf(msr::Symbol,cmp::PropDict)
    if haskey(cmp,:std) return UGF(msr,cmp[:std]) end
    if haskey(cmp,:ntw) ntw, id = cmp[:ntw]; return ntw.usr[id][:ugf] end
    return UGF(msr)
end
function src_ugf(msr::Symbol,src::PropDict)
    if haskey(src,:dep) return UGF(msr,[(1.0)get_unit(msr)],[1.0]) end
    if haskey(src,:std) return UGF(msr,src[:std]) end
    if haskey(src,:ntw) ntw, id = src[:ntw]; return ntw.usr[id][:ugf] end
    return UGF(msr)
end
function set_ugf!(ntw::AbstractNetwork)
    msr = ntw.props[:msr]
    for cmp in cmp(ntw) if !haskey(cmp,:ugf) cmp[:ugf] = cmp_ugf(msr,cmp) end end
    for src in src(ntw) if !haskey(src,:ugf) src[:ugf] = src_ugf(msr,src) end end
end
