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
ld = [(0.603u"km",150u"mm^2"),(0.691u"km",150u"mm^2"),(0.809u"km",150u"mm^2"),
      (0.873u"km",150u"mm^2"),(0.952u"km",150u"mm^2"),(1.095u"km",150u"mm^2"),
      (1.272u"km",150u"mm^2"),(1.464u"km",150u"mm^2"),(2.100u"km",150u"mm^2"),
      (0.603u"km",240u"mm^2"),(0.809u"km",240u"mm^2"),(0.873u"km",240u"mm^2"),
      (0.952u"km",240u"mm^2"),(1.095u"km",240u"mm^2"),(1.272u"km",240u"mm^2"),
      (1.464u"km",240u"mm^2"),(1.678u"km",240u"mm^2"),(0.603u"km",500u"mm^2"),
      (0.691u"km",500u"mm^2"),(1.272u"km",500u"mm^2"),(1.464u"km",500u"mm^2"),
      (1.678u"km",500u"mm^2"),(2.000u"km",500u"mm^2"),(2.100u"km",500u"mm^2"),
      (3.000u"km",500u"mm^2"),(4.000u"km",500u"mm^2"),(4.500u"km",500u"mm^2"),
      (5.000u"km",500u"mm^2"),(5.750u"km",500u"mm^2"),(8.500u"km",500u"mm^2")]

# Listing of the failure rate and repair dictionaries in function of diameter
λ = Dict(150u"mm^2" => 0.00743u"1/yr/km", 
         240u"mm^2" => 0.00743u"1/yr/km",
         500u"mm^2" => 0.00945u"1/yr/km")
μ = Dict(150u"mm^2" => 0.0006944u"1/hr", 
         240u"mm^2" => 0.0006944u"1/hr", 
         500u"mm^2" => 0.0006944u"1/hr")

# Solve the cable function
function solve_cable(nld::Tuple,λ::Dict,μ::Dict)
    nl, nd = nld[1], nld[2]
    std = STD()
    if mss
        add_states!(std, power = [(Inf)u"MW",0.0u"MW",0.0u"MW",0.0u"MW"],
                         init  = [1.0,0.0,0.0,0.0])
        add_transitions!(std, states = [(1,2),(2,1)],
                              rate = [λ[nd]*nl,μ[nd]])
        add_transitions!(std, states = [(1,3),(3,1)], 
                              rate = [0.00168u"1/yr",0.000926u"1/hr"])
        add_transitions!(std, states = [(1,4),(4,1)], 
                              rate = [0.00168u"1/yr",0.000926u"1/hr"])
    else
        λ̄ = λ[nd] * nl + 2 * 0.00168u"1/yr"
        μ̄ = λ̄ / (λ[nd] * nl/μ[nd] + 2 * 0.00168u"1/yr"/ 0.000926u"1/hr")
        add_states!(std, power = [(Inf)u"MW",0.0u"MW"],
                         init  = [1.0,0.0])
        add_transitions!(std, states = [(1,2),(2,1)],
                              rate = [λ̄, μ̄])
    end
    solve!(std, SteadyStateProcess())
    return std
end

# Build a dictionary for the different cables
stdᶜᵇˡ = Dict((nld) => solve_cable(nld,λ,μ) for nld in ld)
