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

# Initializing the state transition diagram of a homogeneous semi markov process to test
# the model

std = STD()
add_states!(std, name  = ["normal operation state","fault operation state f1","Primary post-fault state f1","Secondary post-fault state f1", "Fault operation state for f1", "Primary post-fault state for f2", "secondary post-fault state for f2", "fault operation state for f3", "Primary post-fault state for f3", "secondary post-fault state for f3"],
                        power = [1.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW"],
                        init  = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

add_transitions!(std, distr = [Weibull(10000.0u"hr", 1.25, 0.33), LogNormal(10.0u"hr", 0.08u"hr", 0.8), LogNormal(13.0u"hr", 0.08u"hr", 0.2), LogNormal(100.0u"hr", 0.08u"hr"), Weibull(1000.0u"hr", 10), Weibull(8000.0u"hr", 1.5, 0.33), LogNormal(10.0u"hr", 0.08u"hr", 0.8), LogNormal(13.0u"hr", 0.08u"hr", 0.2), LogNormal(100u"hr", 0.08u"hr"), Weibull(1000.0u"hr", 10), Weibull(12500.0u"hr", 1.75, 0.33), LogNormal(10.0u"hr", 0.08u"hr", 0.8), LogNormal(13.0u"hr", 0.08u"hr", 0.2), LogNormal(100.0u"hr", 0.08u"hr"), Weibull(1000.0u"hr", 10.0)],
                      states = [(1,2),(2,3),(2,4),(4,3),(3,1),(1,5),(5,6),(5,7),(7,6),(6,1),(1,8),(8,9),(8,10),(10,9),(9,1)])
                      
# solve the network
solve!(std, cls , tsim = 450000u"hr", dt = 300u"hr", tol = 1e-8)

Φ1 = _MSS.get_prop(std, 1, :prob)
Φ2 = _MSS.get_prop(std, 2, :prob)
Φ3 = _MSS.get_prop(std, 3, :prob)
Φ4 = _MSS.get_prop(std, 4, :prob)
Φ5 = _MSS.get_prop(std, 5, :prob)
Φ6 = _MSS.get_prop(std, 6, :prob)
Φ7 = _MSS.get_prop(std, 7, :prob)
Φ8 = _MSS.get_prop(std, 8, :prob)
Φ9 = _MSS.get_prop(std, 9, :prob)
Φ10 = _MSS.get_prop(std, 10, :prob)


using Plots

t = 0:300:450000

plot(t, Φ2)
plot(t, Φ3)
plot(t, Φ4)
plot(t, Φ5)
plot(t, Φ6)
plot(t, Φ7)
plot(t, Φ8)
plot(t, Φ9)
plot(t, Φ10)
plot(t, Φ1)
# 
#  