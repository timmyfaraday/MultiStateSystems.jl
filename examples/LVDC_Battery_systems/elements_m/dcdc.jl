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
stdmᵈᶜᵈᶜ = STD() 

# add the states to the std
add_states!(stdmᵈᶜᵈᶜ, name  = ["available", "unavailable_igbt1", "unavailable_igbt2"],
                   power = [10.0u"MW", 0.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0, 0.0])

# add the transitions to the std
add_transitions!(stdmᵈᶜᵈᶜ, states = [(1,2),(2,1),(1,3),(3,1)],
                      distr = [ Exponential(2*(55.6*gamma(1+1/2.31))u"yr"),
                                Exponential((exp(log(4)+(0.3^2)/2))u"d"),
                                Exponential(2*(21.1*gamma(1+1/3.24))u"yr"),
                                Exponential((exp(log(4)+(0.3^2)/2))u"d")])


# solve the std
return stdmᵈᶜᵈᶜ