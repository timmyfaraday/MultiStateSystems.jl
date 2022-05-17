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
stdᵃᶜ = STD()

# add the states to the std
add_states!(stdᵃᶜ, name  = ["available", "unavailable"],
                   power = [1.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdᵃᶜ, rate = [0.0u"1/hr"    3.142e-4u"1/hr"
                                0.8058u"1/hr" 0.0u"1/hr"])

# return the std
return stdᵃᶜ