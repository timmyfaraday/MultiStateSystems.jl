################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a multi-state Markov model considering a subsystem of three
identical and independent channels taken from:

> Reliability of Safety-Critical System - Theory and Applications by Marvin 
  Rausand (2014)

Correspondig slides can be found at:
https://www.ntnu.edu/documents/624876/1277046207/SIS+book+-+chapter+05+-+Introduction+to+Markov+methods/d98c074f-a48e-45b1-a016-668350efbbe7

States:
- 1 : Three channels are functioning
- 2 : Two channels are functioning and one is failed
- 3 : One channel is functioning and two are failed
- 4 : Three channels are failed
"""

# load pkgs
using Unitful
using MultiStateSystems

# initialize the state-transition diagram
std²³ = STD()

# add the states to the std 
add_states!(std²³, init = [1.0, 0.0, 0.0, 0.0])

# add the transitions to the std 
add_transition!(std²³,  states = (1,2),
                        rate   = 3λ)
add_transitions!(std²³, states = [(2,3), (2,1)],
                        rate   = [2λ, μ])
add_transitions!(std²³, states = [(3,4), (3,1)],
                        rate   = [λ, μ])
add_transition!(std²³,  states = (4,1),
                        rate   = μ)

# return the std
return std²³