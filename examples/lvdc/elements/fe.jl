################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Unitful
using MultiStateSystems

# initialize the state-transition diagram corresponding to the front-end ac/dc 
# converter
stdᶠᵉ = STD()

# add the states to the std
add_states!(stdᶠᵉ, name  = ["available", "unavailable"],
                   power = [1.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdᶠᵉ, rate = [0.0u"1/hr"    2.9666667e-6u"1/hr"
                                0.1400u"1/hr" 0.0u"1/hr"])

# return the std
return stdᶠᵉ