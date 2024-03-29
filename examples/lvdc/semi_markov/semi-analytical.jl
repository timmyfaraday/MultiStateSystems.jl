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
using Graphs
using QuadGK

# pkg const
const _MSS = MultiStateSystems
const Graphs  = Graphs

# constants
dt      = 0.1u"hr"
tsim    = 4500.0u"hr"
t       = zero(dt):dt:tsim

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

# # solve the std - Analytical
dτ = dt / 10
dst1 = _MSS.get_prop(std,Graphs.Edge(1,2),:distr)
dst2 = _MSS.get_prop(std,Graphs.Edge(2,3),:distr)

PA1 = 1.0 .- _MSS.cdf.(dst1,t,zero(dt))
PA2 = [quadgk(x -> _MSS.pdf.(dst1,x,zero(nt)) * (1.0 - _MSS.cdf(dst2,nt.-x,nt)),zero(nt),nt,rtol=1e-8)[1] for nt in t]
PA22= [sum(dτ * _MSS.weights(length(zero(nt):dτ:nt))[ni] * _MSS.pdf.(dst1,nl,zero(nt)) * (1.0 - _MSS.cdf(dst2,nt.-nl,nt)) for (ni,nl) in enumerate(zero(nt):dτ:nt)) for nt in t]
PA3 = 1.0 .- (PA1 .+ PA2)



# # solve the std - Semi-Markov
# # setting for a specific analysis
# cls = SemiMarkovProcess()

# solve the std - SEMIMarkov
solve!(std, cls, tsim = tsim, dt = dt)

# # probabilities - SEMIMarkov 
PS1 = _MSS.get_prop(std,1,:prob)
PS2 = _MSS.get_prop(std,2,:prob)
PS3 = _MSS.get_prop(std,3,:prob)


plot(t, [PS1, PS2, PS3],
         label=["State 1" "State 2" "State 3"],
         xlabel="time",
         ylabel="probability")

plot(t, PA1.+PA2,
         label="State 1+ State 2",
         xlabel="time",
         ylabel="probability")
    
ϵ = 1e-16;

plot(t, [PA1.-PS1.+ϵ, PA2.-PS2.+ϵ, PA3.-PS3.+ ϵ],
         label=["State 1" "State 2" "State 3"],
         xlabel="time",
         yscale = :log,
         ylabel="probability difference")

# plot!(t, PA2.-PS2.+ϵ,
#          label="State 2",
#          xlabel="time",
#          yscale = :log,
#          ylabel="probability difference")

# plot!(t, PA3.-PS3.+ ϵ,
#          label="State 3",
#          xlabel="time",
#          yscale = :log,
#          ylabel="probability difference")

#         # plot all probabilities
# plot(t, [PA1, PA2, PA3],
#       label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
#       xlabel="time",
#       ylabel="probability")
        
#          # plot all probabilities
# plot!(t, [PS1, PS2, PS3],
#         label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
#         xlabel="time",
#         ylabel="probability")

# plot(t, h4[2](ustrip(t))-h3[2](ustrip(t)),
#          label=["H array diff for State 2"],
#          xlabel="time",
#          ylabel="rate (1/t)")


# plot the reliability
  plot(_MSS.get_prop(std,:time),
        _MSS.get_prop(std,1,:prob),
        label="reliability",
        xlabel="time",
        ylabel="probability")