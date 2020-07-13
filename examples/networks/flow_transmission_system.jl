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
  Industrial Managers by A. Lisnianski, I. Frenkel, and Y. Ding, pg. 11 (2010)

The oil flow is transmitted from a source S to a user U (Fig. 1).

        n₁  pipe 1  n₂
        |----///----|   pipe 3  n₃
    S⋅--|           |----///----|--⋅ U
        |----///----|
            pipe 2
Fig. 1: Flow transmission system.

The pipes' performance is measured by their transmission capacity (kg/min)
(Table 1). The associated probabilities are given Table 2.

Table 1: Transmission capacity of the pipes in the flow transmission system.
| pipe 1               | pipe 2               | pipe 3               |
|----------------------|----------------------|----------------------|
| g¹₁ = 0    kg/hr     | g²₁ = 0    kg/hr     | g³₁ = 0    kg/hr     |
| g¹₂ = 1500 kg/hr     | g²₂ = 2000 kg/hr     | g³₂ = 1800 kg/hr     |
|                      |                      | g³₃ = 4000 kg/hr     |

Table 2: Associated probabilities of the pipes in the flow transmission system.
| pipe 1               | pipe 2               | pipe 3               |
|----------------------|----------------------|----------------------|
| p¹₁ = 0.2            | p²₁ = 0.4            | p³₁ = 0.1            |
| p¹₂ = 0.8            | p²₂ = 0.6            | p³₂ = 0.2            |
|                      |                      | p³₃ = 0.7            |
"""

# Load Pkgs
using Unitful
using MultiStateSystems

# Initialize the network corresponding to the flow transmission system.
ntw = Network()

# Add the user, source and components to the network ntw.
add_user!(ntw,node = 3)
add_source!(ntw,node = 1)
add_components!(ntw,edge = [(1,2),(1,2),(2,3)],
                    name = ["pipe 1","pipe 2","pipe 3"],
                    std  = [STD(flow = [0u"kg/hr",1500u"kg/hr"],
                                prob = [0.2,0.8]),
                            STD(flow = [0u"kg/hr",2000u"kg/hr"],
                                prob = [0.4,0.6]),
                            STD(flow = [0u"kg/hr",1800u"kg/hr",4000u"kg/hr"],
                                prob = [0.1,0.2,0.7])])

# Solve the problem
@time solve!(ntw)
