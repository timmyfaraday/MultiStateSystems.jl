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
stdⁱⁿᵛ¹ = STD() 

# add the states to the std
add_states!(stdⁱⁿᵛ¹, name  = ["available", "unavailable_weibull"],
                   power = [10.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdⁱⁿᵛ¹, states = [(1,2),(2,1)],
                      distr = [ Weibull(19.2u"yr", 2.64),
                                LogNormal(log(4)u"d", 0.3u"d")])

# solve the std
return stdⁱⁿᵛ¹