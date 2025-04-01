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
stdᵃᶜᵈᶜ = STD() 

# add the states to the std
add_states!(stdᵃᶜᵈᶜ, name  = ["available", "unavailable_weibull"],
                   power = [(Inf)u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

add_transitions!(stdᵃᶜᵈᶜ, states = [(1,2),(2,1)],
            distr = [Weibull(20.0u"yr", 2.38),
                     LogNormal(log(4)u"d", 0.3u"d")])

# solve the std
return stdᵃᶜᵈᶜ