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

  WT₁         WT₂
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

# Include the stochastic models for the wind turbines and cables
include("../elements/wind_turbine.jl")
include("../elements/sea_cable.jl")

# Initialize the network corresponding to the wind farm.
ntw = Network()

# Add the user, sources and components to the network ntw.
add_user!(ntw, node = 8)
add_sources!(ntw, node = 1:6, std = stdʷᵗᵒ, dep = true)
add_components!(ntw, node = 1:6, std = stdʷᵗʳ)
add_components!(ntw, edge = [(1,2),(2,7),(3,4),(4,7),(5,6),(6,7),(7,8)],
                     std = stdᶜᵇˡ)

# Solve the problem
solve!(ntw, type = :steady)