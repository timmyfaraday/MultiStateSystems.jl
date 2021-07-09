################################################################################
#  Copyright 2020, Tom Van Acker                                               #
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
stdᵃᵈ = STD()

# add the states to the std
add_states!(stdᵃᵈ, name  = ["available", "unavailable"],
                   power = [1.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdᵃᵈ, rate = [0.0u"1/hr"    12.8928e-6u"1/hr"
                                0.0833u"1/hr" 0.0u"1/hr"])

# return the std
return stdᵃᵈ