################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
# Measure (abbr: msr)
"""
# sets
const MsrSet   = Set{Symbol}([:flow,:power])
const UnitDict = Dict{Symbol,_UF.FreeUnits}(:flow   => u"m^3/hr",
                                            :power  => u"MW")
const MaxVal   = Dict{Symbol,_UF.Quantity}( :flow   => (Inf)u"m^3/hr",
                                            :power  => (Inf)u"MW")
# functions
get_max(msr::Symbol) = MaxVal[msr]
get_unit(msr::Symbol) = UnitDict[msr]
check_unit(msr::Symbol,unit::_UF.FreeUnits) =
    _UF.dimension(UnitDict[msr]) == _UF.dimension(unit)

get_msr(std::AbstractSTD) = haskey(std.props,msr) ? std.props[:msr] : Set() ;
get_msr(ntw::AbstractNetwork) = haskey(ntw.props,msr) ? ntw.props[:msr] : Set() ;
set_msr(ntw::AbstractNetwork) =
    for ne in elements(ntw)
        if haskey(kwargs,:std) push!(get_msr(ntw),get_msr(ne[:std])...) end
        if haskey(kwargs,:ntw) push!(get_msr(ntw),get_msr(ne[:ntw][1])...) end
    end
