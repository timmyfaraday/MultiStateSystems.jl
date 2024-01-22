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

# initialize the state-transition diagram corresponding to the dc/dc converter capacitance
stdmᶜᵃᵖᵈᶜ = STD() 

# add the states to the std
add_states!(stdmᶜᵃᵖᵈᶜ, name  = ["available", "unavailable"],
                   power = [10.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])


# add the transitions to the std
add_transitions!(stdmᶜᵃᵖᵈᶜ, states = [(1,2),(2,1)],
                      distr = [ Exponential((31.65*gamma(1+1/8.5))u"yr"), 
                                Exponential((exp(log(4)+(0.3^2)/2))u"d")])

# solve the std
return stdmᶜᵃᵖᵈᶜ