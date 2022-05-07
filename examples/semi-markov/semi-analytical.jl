################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a semi-Markov model taken from:

> Mathematical formulation and numerical treatment based on transition frequency 
  densities and quadrature methods for non-homogeneous semi-Markov processes by 
  Moura and Droguett (2009)

States:
- 1 : normal operation state
- 2 : degradation state
- 3 : failed state
"""

# load pkgs
using Plots
using Unitful
using UnitfulRecipes
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SemiMarkovProcess()

# initialize the state-transition diagram
std = STD()

# add the states to the std 
add_states!(std, name  = ["normal operation state","degradation state","failed state"],
                 init  = [1.0,0.0,0.0])

# add the transitions to the std 
add_transitions!(std, distr = [Exponential(1000.0u"hr"), Weibull(250.0u"hr", 1.5)],
                      states = [(1,2),(2,3)])
                      
# solve the std
solve!(std, cls, tsim = 4500u"hr", dt = 3u"hr")

# plot all probabilities
plot(_MSS.get_prop(std, :time), 
        [_MSS.get_prop(std, ns, :prob) for ns in _MSS.states(std)],
        label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
        xlabel="time",
        ylabel="probability")

# plot the reliability
plot(_MSS.get_prop(std, :time),
        _MSS.get_prop(std, 1, :prob) + _MSS.get_prop(std, 2, :prob),
        label="reliability",
        xlabel="time",
        ylabel="probability")