################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a multi-state Markov model for an optical monitoring system
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
add_states!(std, name  = ["normal operation state","reinstallation state","restore state","halted state"],
                 init  = [1.0,0.0,0.0,0.0])

# add the transitions to the std 
add_transitions!(std, states = [(1,1),(1,2),(2,1),(2,3),(3,1),(3,4)],
                      distr = [ Exponential(10000.0u"hr", t -> -0.00004u"1/hr" * t + 0.6), 
                                Weibull(150.0u"hr", 1.36, t -> 0.00004u"1/hr" * t + 0.4), 
                                Exponential(20.0u"hr", 0.82), 
                                LogNormal(2.5u"hr", 0.25u"hr", 0.18), 
                                Exponential(20.0u"hr", 0.62), 
                                LogNormal(4.0u"hr", 0.4u"hr", 0.38)])

# return the std
return std