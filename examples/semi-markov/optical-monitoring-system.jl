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

# Initializing the state transition diagram of the application example as can be
# found in "Mathematical formulation and numerical treatment based on transition
# frequency densities and quadrature methods for non-homogeneous semi-Markov processes"
# by Moura and Droguett.

std = STD()
add_states!(std, name  = ["normal operation state","reïnstallation state","restore state","halted state"],
                        power = [1.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW"],
                        init  = [1.0,0.0,0.0,0.0])

add_transitions!(std, distr = [Exponential(10000.0u"hr", x -> -0.00004x + 0.6), Weibull(150.0u"hr", 1.36, x -> 0.00004x + 0.4), Exponential(20u"hr", 0.82), LogNormal(2.5u"hr", 0.25u"hr", 0.18), Exponential(20u"hr", 0.62), LogNormal(4.0u"hr", 0.4u"hr", 0.38)],
                      states = [(1,1),(1,2),(2,1),(2,3),(3,1),(3,4)])
                      
# solve the network
solve!(std, cls , tsim = 8760u"hr", dt = 6u"hr", tol = 1e-8)

Φ1 = _MSS.get_prop(std, 1, :prob)
Φ2 = _MSS.get_prop(std, 2, :prob)
Φ3 = _MSS.get_prop(std, 3, :prob)
Φ4 = _MSS.get_prop(std, 4, :prob)

using Plots

t = 0:6:8760

plot(t, Φ2)
plot(t, Φ3)
plot(t, Φ4)
plot(t, Φ1)
