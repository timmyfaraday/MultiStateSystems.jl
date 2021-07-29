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
has_msr(kwargs::Iterators.Pairs) =  
    !isempty(collect(intersect(MsrSet,keys(kwargs))))
get_msr(ugf::AbstractUGF) = [ugf.msr]
get_msr(kwargs::Iterators.Pairs) = collect(intersect(MsrSet,keys(kwargs)))[1]
get_msr(std::AbstractSTD) = haskey(std.props,:msr) ? [std.props[:msr]] : [] ;
get_msr(ntw::AbstractNetwork) = haskey(ntw.props,:msr) ? [ntw.props[:msr]] : [] ;
function set_msr!(ntw::AbstractNetwork)
    msr = []
    for ne in elements(ntw)
        if haskey(ne,:ugf) union!(msr,get_msr(ne[:ugf])) end
        if haskey(ne,:std) union!(msr,get_msr(ne[:std])) end
        if haskey(ne,:ntw) union!(msr,get_msr(ne[:ntw][1])) end
    end
    ntw.props[:msr] = msr[1]                                                    # currently, only one msr is supported.
end
