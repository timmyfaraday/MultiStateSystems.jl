################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a flow transmission system with three pipes taken from:

> Multi-State System Reliability Analysis and Optimization for Engineers and
  Industrial Managers by A. Lisnianski, I. Frenkel, and Y. Ding, pg. 11/209 
  (2010)

The oil flow is transmitted from a source S to a user U (Fig. 1).

        n₁  pipe 1  n₂
        |----///----|   pipe 3  n₃
    S⋅--|           |----///----|--⋅ U
        |----///----|
            pipe 2
Fig. 1: Flow transmission system.

The pipes' performance is measured by their transmission capacity (m³/hr)
(Table 1).

Table 1: Transmission capacity of the pipes in the flow transmission system.
| pipe 1               | pipe 2               | pipe 3               |
|----------------------|----------------------|----------------------|
| g¹₁ = 0    m³/hr     | g²₁ = 0    m³/hr     | g³₁ = 0    m³/hr     |
| g¹₂ = 1500 m³/hr     | g²₂ = 2000 m³/hr     | g³₂ = 1800 m³/hr     |
|                      |                      | g³₃ = 4000 m³/hr     |
"""

# load pkgs
using Unitful
using MultiStateSystems

# initialize the state-transition diagrams corresponding with the pipes
stdᵖ¹, stdᵖ², stdᵖ³ = STD(), STD(), STD()

# add the states and transitions to the respective state-transition diagrams
## pipe 1
add_states!(stdᵖ¹, flow = [0u"m^3/hr", 1500u"m^3/hr"],
                   init = [0.0, 1.0])
add_transitions!(stdᵖ¹, states = [(1,2), (2,1)],
                        rate   = [μ¹, λ¹])
## pipe 2
add_states!(stdᵖ², flow = [0u"m^3/hr", 2000u"m^3/hr"],
                   init = [0.0, 1.0])
add_transitions!(stdᵖ², states = [(1,2), (2,1)],
                        rate   = [μ², λ²])
## pipe 3
add_states!(stdᵖ³, flow = [0u"m^3/hr", 1800u"m^3/hr", 4000u"m^3/hr"],
                   init = [0.0, 0.0, 1.0])
add_transitions!(stdᵖ³, states = [(1,2), (2,1), (2,3), (3,2)],
                        rate   = [μ³ᵃ, λ³ᵃ, μ³ᵇ, λ³ᵇ])

# solve the state-transition diagrams 
solve!(stdᵖ¹, SteadyStateProcess())
solve!(stdᵖ², SteadyStateProcess())
solve!(stdᵖ³, SteadyStateProcess())

# initialize the network corresponding to the flow transmission system
ntwᶠᵗˢ = Network()

# add the user, source and components to the network ntwᶠᵗˢ
add_user!(ntwᶠᵗˢ, node = 3)
add_source!(ntwᶠᵗˢ, node = 1)
add_components!(ntwᶠᵗˢ, edge = [(1,2),(1,2),(2,3)],
                        name = ["pipe 1", "pipe 2", "pipe 3"],
                        std  = [stdᵖ¹, stdᵖ², stdᵖ³])

# return the network
return ntwᶠᵗˢ
