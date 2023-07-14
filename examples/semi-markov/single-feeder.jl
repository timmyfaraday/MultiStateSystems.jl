################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Plots
using Unitful
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SemiMarkovProcess()

# Initializing the state transition diagram of a homogeneous semi markov process to test
# the model

std = STD()
add_states!(std, name  = ["normal operation state","fault operation state f1","Primary post-fault state","Secondary post-fault state", "Fault operation state for f2", "Fault operation state for f3"],
                        power = [1.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW"],
                        init  = [1.0,0.0,0.0,0.0,0.0,0.0])

add_transitions!(std, states = [(1,2),(1,5),(1,6),(2,3),(5,3),(6,3),(2,4),(5,4),(6,4),(3,4),(4,1)],
                      distr = [Exponential(4000.0u"hr", 0.3),
                               Exponential(10000.0u"hr", 0.4),
                               Exponential(25000.0u"hr", 0.3),
                               LogNormal(10.0u"hr", 0.08u"hr", 0.8),
                               LogNormal(10.0u"hr", 0.08u"hr", 0.8),
                               LogNormal(10.0u"hr", 0.08u"hr", 0.8),
                               LogNormal(10.0u"hr", 0.08u"hr", 0.2),
                               LogNormal(10.0u"hr", 0.08u"hr", 0.2),
                               LogNormal(10.0u"hr", 0.08u"hr", 0.2),
                               LogNormal(100u"hr", 0.08u"hr"),
                               Weibull(1000.0u"hr", 10)])
                      
                      
# solve the network
solve!(std, cls , tsim = 10000.0u"hr", dt = 10u"hr")

# plot the probabilities
plot(_MSS.get_prop(std, :time), 
        [_MSS.get_prop(std, ns, :prob) for ns in _MSS.states(std)],
        label=reshape([_MSS.get_prop(std, ns, :name) for ns in _MSS.states(std)],1,:),
        xlabel="time",
        ylabel="probability")

# plot the reliability
plot(_MSS.get_prop(std, :time),
        _MSS.get_prop(std, 1, :prob),
        label="semi-Markov",
        xlabel="time",
        ylabel="Availability")


tot_prob = []
ns = 11
tot_prob .+= _MSS.get_prop(std, ns, :prob)
for ns in _MSS.states(std)
        tot_prob .+= _MSS.get_prop(std, ns, :prob)
end
