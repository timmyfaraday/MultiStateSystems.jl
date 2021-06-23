################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a tree-like wind farm with four wind turbines.

The goal of this example is to show the effectiveness of solving tree-like
systems iteratively compared to solving it as a single system. To this end, a
wind farm with four wind turbines is considered (Fig. 1).

  WT₁         WT₂          n₅
   |----///----|----///----|
   n₁    c₁     n₂    c₂     |
                           |--⋅ PCC
  WT₃         WT₄           |
   |----///----|----///----|
   n₃    c₃     n₄    c₄
Fig. 1: Tree-like wind farm with four wind turbines.
"""
# Load Pkgs
using Unitful
using MultiStateSystems

# Initialize the network corresponding to the wind farm.
ntw = Network()

# Initialize the necessary state-transition diagrams
stdʷᵗ = STD(prob  = [0.3, 0.7],
            power = [0.0u"MW", 2.0u"MW"])
stdᶜ  = STD(prob  = [0.1, 0.9],
            power = [0.0u"MW", 3.0u"MW"])

# Add the user, sources and components to the network ntw.
add_user!(ntw, node = 5)
add_sources!(ntw, node = 1:4, std = stdʷᵗ, dep = true)
add_components!(ntw, edge = [(1,2),(2,5),(3,4),(4,5)], std = stdᶜ)

# Solve the problem
solve!(ntw, type = :steady)
