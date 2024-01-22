################################################################################
#  Copyright 2023, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
This file demonstrates the resilience of LVDC systems with a battery, 
against AC side disturbances and outages in function of battery reserve time.
A semi-Markov process strategy is applied to the problem and the state transition
diagram contains three states.

States:
- 1 : Normal operation state
- 2 : Failed AC state, battery is discharging
- 3 : Shutdown state, battery is empty
"""

# load pkgs
using Plots
using Unitful
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SemiMarkovProcess()

bat = 120.0
AC_fail = 6.6u"yr"

function battery_dist(T::Quantity)
# This function determines the transition distribution from an empty battery state to 
# a normal operation state.
        if T == 120.0u"minute"
                return Exponential(51.21u"minute")
        elseif T==60.0u"minute"
                return Exponential(55.4u"minute")
        elseif T==80.0u"minute"
                return Exponential(54.39u"minute")
        elseif T==100.0u"minute"
                return Exponential(53.14u"minute")
        elseif T == 140.0u"minute"
                return Exponential(47.38u"minute")
        end
end     


# initialize the state-transition diagram
stdᵈᶜ = STD()
stdᵃᶜ = STD()

# add the states to the std 
add_states!(stdᵃᶜ, name  = ["normal operation state","failed AC state"],
                init  = [1.0,0.0])

add_states!(stdᵈᶜ, name  = ["normal operation state","failed AC state","Battery depletion state"],
                init  = [1.0,0.0,0.0])
# add the transitions to the std 
add_transitions!(stdᵃᶜ, states = [(1,2),(2,1)],
                distr = [ Exponential(AC_fail), 
                                Weibull(95.196u"minute", 2.1)])

add_transitions!(stdᵈᶜ, states = [(1,2),(2,1),(2,3),(3,1)],
                        distr = [ Exponential(AC_fail), 
                                Weibull(95.196u"minute", 2.1, 0.80),  
                                LogNormal(log(bat/60)u"hr", 0.05u"hr", 0.20),
                                battery_dist((bat)u"minute")])

# solve the stds
solve!(stdᵈᶜ, cls, tsim = 7.0u"hr", dt = 1u"minute", tol=1e-10)
solve!(stdᵃᶜ, cls, tsim = 7.0u"hr", dt = 1u"minute", tol=1e-10)