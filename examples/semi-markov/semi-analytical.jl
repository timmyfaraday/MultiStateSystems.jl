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
using LightGraphs
using QuadGK

# pkg const
const _MSS = MultiStateSystems
const _LG  = LightGraphs

# constants
dt      = 3.0u"hr"
tsim    = 4500.0u"hr"
t       = zero(dt):dt:tsim

# setting for a specific analysis
cls = MarkovProcess()

# initialize the state-transition diagram
std = STD()

# add the states to the std 
add_states!(std, name  = ["normal operation state","degradation state","failed state"],
                 init  = [1.0,0.0,0.0])

# add the transitions to the std 
add_transitions!(std, distr = [Exponential(1000.0u"hr"), Weibull(250.0u"hr", 1.5)],
                      states = [(1,2),(2,3)])

# trapping
_MSS.set_info!(std, 3, :trapping, true)

# # solve the std - Markov
# solve!(std, cls, tsim = tsim, dt = dt)

# # probabilities - Markov 
# PM1 = _MSS.get_prop(std,1,:prob)
# PM2 = _MSS.get_prop(std,2,:prob)
# PM3 = _MSS.get_prop(std,3,:prob)

# # plot all probabilities
# plot(_MSS.get_prop(std, :time), 
#         [PM1, PM2, PM3],
#         label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
#         xlabel="time",
#         ylabel="probability")

# solve the std - Analytical
dτ = dt / 100
dst1 = _MSS.get_prop(std,_LG.Edge(1,2),:distr)
dst2 = _MSS.get_prop(std,_LG.Edge(2,3),:distr)

PA1 = 1.0 .- _MSS.cdf.(dst1,t,t)
PA2 = [quadgk(x -> _MSS.pdf.(dst1,x,x) * (1.0 - _MSS.cdf(dst2,nt.-x,nt)),zero(nt),nt,rtol=1e-8)[1] for nt in t]
PA3 = 1.0 .- (PA1 .+ PA2)

# plot all probabilities
plot(t, [PA1, PA2, PA3],
        label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
        xlabel="time",
        ylabel="probability")

# solve the std - Semi-Markov
# setting for a specific analysis
cls = SemiMarkovProcess()

# solve the std - SEMIMarkov
solve!(std, cls, tsim = tsim, dt = dt)

# probabilities - SEMIMarkov 
PS1 = _MSS.get_prop(std,1,:prob)
PS2 = _MSS.get_prop(std,2,:prob)
PS3 = _MSS.get_prop(std,3,:prob)

# plot all probabilities
plot!(t, [PS1, PS2, PS3],
         label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
         xlabel="time",
         ylabel="probability")







# # plot the reliability
# plot(_MSS.get_prop(std, :time),
#         _MSS.get_prop(std, 1, :prob) + _MSS.get_prop(std, 2, :prob),
#         label="reliability",
#         xlabel="time",
#         ylabel="probability")