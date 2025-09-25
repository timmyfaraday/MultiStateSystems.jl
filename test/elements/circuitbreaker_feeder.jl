################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Unitful
using MultiStateSystems

# pkg constants
const _MSS = MultiStateSystems

# initialize the state-transition diagram corresponding to the front-end ac/dc capacitance
std¹ = STD() 

# add the states to the std
add_states!(std¹, name  = ["A", "U", "V"],
    power = [(Inf)u"MW", 10.0u"MW", 0.0u"MW"],
    init  = [1.0, 0.0, 0.0])

# add the transitions to the std
add_transitions!(std¹, states = [(1,2),(1,3),(2,1),(3,1)],
    distr = [   Weibull(15.0u"yr", 3.81, 0.9),
                Weibull(15.0u"yr", 3.81, 0.1), 
                LogNormal(log(10.0)u"d", 0.2u"d"),
                LogNormal(log(10.0)u"d", 0.2u"d")])

# initialize the state-transition diagram corresponding to the front-end ac/dc capacitance
std² = STD() 

# add the states to the std
add_states!(std², name  = ["A", "U1", "V", "U2"],
    power = [(Inf)u"MW", 10.0u"MW", 0.0u"MW", 10.0u"MW"],
    init  = [1.0, 0.0, 0.0, 0.0])

# add the transitions to the std
add_transitions!(std², states = [(1,2),(1,3),(2,1),(3,4),(4,1)],
    distr = [   Weibull(15.0u"yr", 3.81, 0.9),
                Weibull(15.0u"yr", 3.81, 0.1), 
                LogNormal(log(10.0)u"d", 0.2u"d"),
                LogNormal(log(2.0)u"hr", 0.25u"hr"),
                LogNormal(log(9.914)u"d", 0.2u"d")])

# return the std
return std¹, std²