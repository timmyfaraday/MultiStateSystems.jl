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
stdᵈᶜᵈᶜ = STD() 

# add the states to the std
add_states!(stdᵈᶜᵈᶜ, name  = ["available", "unavailable_igbt1", "unavailable_igbt2"],
                   power = [10.0u"MW", 0.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0, 0.0])

# add the transitions to the std
add_transitions!(stdᵈᶜᵈᶜ, states = [(1,2),(2,1),(1,3),(3,1)],
                      distr = [ Weibull(55.6u"yr", 2.31, t -> 0.00004u"1/d" * t + 0.5),
                                LogNormal(log(4)u"d", 0.3u"d"),
                                Weibull(21.1u"yr", 3.24, t -> -0.00004u"1/d" * t + 0.5),
                                LogNormal(log(4)u"d", 0.3u"d")])


# solve the std
return stdᵈᶜᵈᶜ