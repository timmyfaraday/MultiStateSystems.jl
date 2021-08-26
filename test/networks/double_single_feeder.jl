################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a double single feeder section taken from:

> Industrial Energy System Availability Management by T. Van Acker, pg. 92 (2020)

Power is transmitted from sources S₁/S₂ to users U₁/U₂ using two independent 
single feeders (Fig. 1)

      n₁     fdr 1    n₃
    S₁ ⋅-□----///----□-⋅ U₁

    S₂ ⋅-□----///----□-⋅ U₂
      n₂     fdr 2    n₄
Fig. 1: Double single feeder section
"""

# load pkgs
using Unitful
using MultiStateSystems

# initialize the network corresponding with the double single feeder section
ntwᵈˢᶠ = Network()

# add the user, source and components to the network ntwᵈˢᶠ
add_users!(ntwᵈˢᶠ, node = [3,4])
add_sources!(ntwᵈˢᶠ, node = [1,2])
add_components!(ntwᵈˢᶠ, edge = [(1,3),(2,4)],
                        name = ["fdr 1", "fdr 2"],
                        std  = [solvedSTD(power = [0.0u"MW", 1.0u"MW"],
                                          prob  = [0.1,0.9]),
                                solvedSTD(power = [0.0u"MW", 1.0u"MW"],
                                          prob  = [0.1,0.9])])

# return the network
return ntwᵈˢᶠ