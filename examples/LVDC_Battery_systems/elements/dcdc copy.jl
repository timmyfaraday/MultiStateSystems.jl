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

# initialize the state-transition diagram corresponding to the front-end ac/dc converter
std = STD() 

# add the states to the std
add_states!(std, name  = ["available", "unavailable_igbt1", "unavailable_igbt2"],
                   power = [10.0u"MW", 0.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0, 0.0])

# add the transitions to the std
add_transitions!(std, states = [(1,2),(2,1),(1,3),(3,1)],
                      distr = [ Exponential(2*55.6u"yr"),
                                Exponential(10.0u"d"),
                                Exponential(2*21.1u"yr"),
                                Exponential(10.0u"d")])


# solve the std
return std