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
- 2 : reinstallation state
- 3 : restoration state
- 4 : halted state
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
                      
# solve the std
solve!(std, cls, tsim = 8760.0u"hr", dt = 5.84u"hr")

# plot the probabilities
plot(_MSS.get_prop(std, :time), 
        [_MSS.get_prop(std, ns, :prob) for ns in _MSS.states(std)],
        label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
        xlabel="time",
        ylabel="probability")

# plot the reliability
plot!(_MSS.get_prop(std, :time),
        _MSS.get_prop(std, 1, :prob),
        label="reliability",
        xlabel="time",
        ylabel="probability")