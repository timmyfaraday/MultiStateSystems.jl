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
stdᵇᵃᵗ = STD()

# add the states to the std
add_states!(stdᵇᵃᵗ, name  = ["available", "unavailable"],
                  power = [10.0u"MW", 0.0u"MW"],
                  init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdᵇᵃᵗ,  states = [(1,2),(2,1)],
                        distr = [Exponential(55.0u"yr"),
                                LogNormal(log(4)u"d", 0.3u"d")])
# return the std
return stdᵇᵃᵗ