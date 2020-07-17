################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# Load Pkgs
using Unitful
using MultiStateSystems

# Listing of possible combination of lengths and diameters
ld = [(10u"km",240u"mm^2"),(20u"km",240u"mm^2")]

# Listing of the failure rate and repair dictionaries in function of diameter
λ = Dict(240u"mm^2" => 0.000795u"1/yr/km")
μ = Dict(240u"mm^2" => 0.075u"1/hr")

# Solve the cable function
function solve_cable(nld::Tuple,λ::Dict,μ::Dict)
    nl, nd = nld[1], nld[2]
    std = STD()
    add_states!(std, flow = [(Inf)u"MW",0.0u"MW",0.0u"MW",0.0u"MW"],
                     init = [1.0,0.0,0.0,0.0])
    add_transitions!(std, states = [(1,2),(2,1)],
                          rate = [λ[nd]*nl,μ[nd]])
    add_transitions!(std, states = [(1,3),(3,1)],                   #cable termination1
                          rate = [0.00168u"1/yr",0.000926u"1/hr"])
    add_transitions!(std, states = [(1,4),(4,1)],                   #cable termination2
                          rate = [0.00168u"1/yr",0.000926u"1/hr"])
    solve!(std, 1.0u"yr")
    return std
end

# Build a dictionary for the different cables
stdᶜᵇˡ = Dict((nld) => solve_cable(nld,λ,μ) for nld in ld)
