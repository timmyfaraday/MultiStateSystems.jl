################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a multi-state Markov model for a optical monitoring system
taken from:

> A Multi-State Markov Model for a Short-Term Reliability Analysis of a Power
  Generating Unit by A. Lisnianski, D. Elmakias, and H. Ben Haim (2012)

"""

# load pkgs
using Unitful
using MultiStateSystems



# initialize the state-transition diagram
std = STD()

# add the states to the std
add_states!(std, name  = ["available", "unavailable"],
                   power = [1.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])

# add the transitions to the std
add_transitions!(stdᵃᵈ, rate = [0.0u"1/hr"    12.8928e-6u"1/hr"
                                0.0833u"1/hr" 0.0u"1/hr"])

# return the std
return stdᵃᵈ