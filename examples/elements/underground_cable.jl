################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# Load Pkgs
using Unitful
using MultiStateSystems

# Initialize the state-transition diagram corresponding to the transformer
stdᵗʳ = STD()

# Add the states to the std
add_states!(stdᵗʳ, init  = [1.0, 0.0, 0.0, 0.0],
                   φinit = [5.0u"yr", 0.0, 0.0, 0.0],
                   power = [20.0u"MW", 0.0u"MW", 0.0u"MW", 0.0u"MW"])

# Add the transitions to the std
add_transitions!(stdᵗʳ, states = [(1,2), (2,3), (2,4), (4,3), (3,1)],
                        type   = [:f, :r, :r, :r, :cpm],
                        distr  = [𝑾(20.0u"yr",4.0),
                                  𝑫(0.1u"s",0.99),
                                  𝑫(0.3u"s",0.01),
                                  𝑾(5.0u"hr",3.0),
                                  𝑾(2.0u"d",3.0)])

# Solve the STD
# solve!(std, 1.0u"yr", alg = :vanacker_process)