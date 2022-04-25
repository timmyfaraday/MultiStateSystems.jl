################################################################################
#  Copyright 2020, Tom Van Acker                                               #
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

std = STD()
add_states!(std, name  = ["normal operation state","degradation state","failed state"],
                        power = [1.0u"MW", 1.0u"MW", 0.0u"MW"],
                        init  = [1.0,0.0,0.0])

add_transitions!(std, distr = [Exponential(1000u"hr"), Weibull(250u"hr", 1.5)],
                      states = [(1,2),(2,3)])
                      
# solve the network
solve!(std, cls)

println("test")