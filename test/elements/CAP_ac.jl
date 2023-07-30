################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Unitful
using MultiStateSystems
using Revise

const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SemiMarkovProcess()

# initialize the state-transition diagram corresponding to the front-end ac/dc capacitance
stdᶜᵃᵖᵃᶜ = STD() 

# add the states to the std
add_states!(stdᶜᵃᵖᵃᶜ, name  = ["available", "unavailable"],
                   power = [10.0u"MW", 0.0u"MW"],
                   init  = [1.0, 0.0])


# add the transitions to the std
add_transitions!(stdᶜᵃᵖᵃᶜ, states = [(1,2),(2,1)],
                      distr = [ Exponential(5.0u"yr"),
                                Exponential(48.0u"d")])
                              
                                # LogNormal(3.0u"yr", 1.0u"yr")

# dt = 4u"hr"
# t = 0u"hr":dt:8760.0u"hr";

# @time U = set_U(stdᶜᵃᵖᵃᶜ, t, 1e-8)

_MSS.solveP!(stdᶜᵃᵖᵃᶜ, cls, tsim = 20.0u"hr", dt = 4u"hr", tol = 1e-9);

t = 0.0u"hr":4u"hr":1.0u"yr";



# solve the std
return stdᶜᵃᵖᵃᶜ