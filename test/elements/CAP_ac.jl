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

# pkg constants
const _MSS = MultiStateSystems

# initialize the state-transition diagram corresponding to the front-end ac/dc capacitance
std = STD() 

# add the states to the std
add_states!(std, name  = ["available", "unavailable"],
                 power = [10.0u"MW", 0.0u"MW"],
                 init  = [1.0, 0.0])


# add the transitions to the std
add_transitions!(std, states = [(1,2),(2,1)],
                      distr  = [Exponential(5.0u"yr"),
                                Exponential(48.0u"d")])

# solve the std
return std