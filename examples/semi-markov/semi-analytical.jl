################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Unitful
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SemiMarkovProcess()

# Initializing the state transition diagram of the semi analytical example as can be
# found in "Mathematical formulation and numerical treatment based on transition
# frequency densities and quadrature methods for non-homogeneous semi-Markov processes"
# by Moura and Droguett.

std = STD()
add_states!(std, name  = ["normal operation state","degradation state","failed state"],
                        power = [1.0u"MW", 1.0u"MW", 0.0u"MW"],
                        init  = [1.0,0.0,0.0])

add_transitions!(std, distr = [Exponential(1000.0u"hr"), Weibull(250.0u"hr", 1.5)],
                      states = [(1,2),(2,3)])
                      
# solve the network
solve!(std, cls , tsim = 4500u"hr", dt = 3u"hr", tol = 1e-8)

Φ1 = _MSS.get_prop(std, 1, :prob);
Φ2 = _MSS.get_prop(std, 2, :prob);
Φ  = _MSS.get_prop(std, 1, :prob) + _MSS.get_prop(std, 2, :prob);

using Plots

t = 0:3:4500

plot(t, Φ1)
plot(t, Φ2)
plot(t, Φ)
