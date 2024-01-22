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
using SpecialFunctions

# initialize the state-transition diagram corresponding to the front-end ac/dc converter
stdmᵃᶜᵈᶜ = STD() 

# add the states to the std
add_states!(stdmᵃᶜᵈᶜ, name  = ["available", "unavailable"],
                   power = [10.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

add_transitions!(stdmᵃᶜᵈᶜ, states = [(1,2),(2,1)],
            distr = [   Exponential((20.0*gamma(1+1/2.38))u"yr"),
                        Exponential((exp(log(4)+(0.3^2)/2))u"d")])

# solve the std
return stdmᵃᶜᵈᶜ