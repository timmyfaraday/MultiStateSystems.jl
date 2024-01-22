################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a semi-analytical example of a multi-state semi-Markov model
taken from:

> Mathematical formulation and numerical treatment based on transition frequency
densities and quadrature methods for non-homogeneous semi-Markov processes by 
M. das Chagas Moura, and E. L. Droguett (2009)
"""

# load pkgs
using Unitful
using MultiStateSystems

# initialize the state-transition diagram
std = STD()

# add the states to the std
add_states!(std, name  = ["operational", "non-critical degradation", "failed"],
                 init  = [1.0, 0.0, 0.0])

# add the transitions to the std
add_transitions!(std, states = [(1,2),(2,3)],
                      distr  = [Exponential(1000.0u"hr"),
                                Weibull(250.0u"hr",1.5)])

# return the std
return std