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

A UGF constructor for a specific measure `msr` based on a given probability
vector `prb` and value vector `val`.

This function automatically reduces the state-space to where only unique values
and associated probabilites remain.

# Example
```julia-repl
julia> ugfᵍᵉⁿ = UGF(:flow, [0.0u"MW",0.0u"MW",2.0u"MW"], [0.1,0.2,0.7])
julia> isequal(ugfᵍᵉⁿ,[0.0u"MW",2.0u"MW"])
true
```
"""
struct UGF <: AbstractUGF
    msr::Symbol
    val::Vector
    prb::Vector

    UGF(msr,val,prb) = new(msr,reduce(val,prb)...)
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
julia> stdᵍᵉⁿ = STD(prob = [0.1,0.2,0.7],
                    flow = [0.0u"MW",0.0u"MW",2.0u"MW"])
julia> ugfᵍᵉⁿ = UGF(:flow, stdᵍᵉⁿ)
julia> isequal(ugfᵍᵉⁿ,[0.0u"MW",2.0u"MW"])
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
    for cmp in cmp(ntw) cmp[:ugf] = cmp_ugf(msr,cmp) end
    for src in src(ntw) src[:ugf] = src_ugf(msr,src) end
end
