################################################################################
#  Copyright 2022, Tom Van Acke, Glenn Emmers                                  #
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
stdᵇ = STD()

# add the states to the std
add_states!(stdᵇ, name  = ["available", "unavailable"],
                  power = [1.0u"MW", 0.0u"MW"],
                  init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdᵇ, rate = [0.0u"1/hr"    2.07222222e-6u"1/hr"
                               0.0658u"1/hr" 0.0u"1/hr"])

# return the std
return stdᵇ