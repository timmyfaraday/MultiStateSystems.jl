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
wind farm with six wind turbines is considered (Fig. 1).

  WT₁         WT₂           n₃
   |----///----|----///----|
   n₁    c₁     n₂    c₂     |
                           |--⋅ PCC
  WT₃         WT₄           |
   |----///----|----///----|
   n₅    c₃     n₄    c₄
Fig. 1: Tree-like wind farm with four wind turbines.
"""

# load Pkgs
using Unitful
using MultiStateSystems

# include the state-transition diagrams for the wind turbines and cables
stdʷᵗᵒ = STD(prob = [0.3,0.7], power = [0.0u"MW",2.0u"MW"])
stdᶜᵇˡ = STD(prob = [0.1,0.9], power = [0.0u"MW",4.0u"MW"])

# initialize the network corresponding to the wind farm.
ntwʷᶠ = Network()

# add the user, sources and components to the network ntw.
add_user!(ntwʷᶠ, node = 3)
add_sources!(ntwʷᶠ, node = [1,2,4,5], std = stdʷᵗᵒ, dep = true)
add_components!(ntwʷᶠ, edge = [(1,2),(2,3),(3,4),(4,5)],
                       std  = stdᶜᵇˡ)

# return the network
return ntwʷᶠ