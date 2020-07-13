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

# Initialize the state-transition diagram corresponding to the sub-sea cable.
stdᶜᵇˡ = STD()

# Add the states to the std
add_states!(stdᶜᵇˡ, flow = [0.0u"MW",(Inf)u"MW"],
                    init = [0.0,1.0])

# Add the transitions to the std
add_transitions!(stdᶜᵇˡ, rate = [0.000u"1/yr" 6.083u"1/yr"
                                 0.075u"1/yr" 0.000u"1/yr"])

# Solve the problem
solve!(stdᶜᵇˡ, 1.0u"yr")
