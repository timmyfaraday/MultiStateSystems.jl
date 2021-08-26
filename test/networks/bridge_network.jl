################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
Example considering a bridge network.
"""

# load pkgs 
using Unitful
using MultiStateSystems

# include the state-transition diagrams for the pipes
stdᵖ = solvedSTD(prob = [0.1, 0.9], flow = [0.0u"m^3/hr", 1.0u"m^3/hr"])

# initialize the network corresponding with the bridge network
ntwᵇⁿ = Network()

# add the sources, components and users to the network 
add_source!(ntwᵇⁿ, node = 1)
add_components!(ntwᵇⁿ, edge = [(1,2), (2,4), (1,3), (3,4)],
                       std  = stdᵖ)
add_bidirectional_component!(ntwᵇⁿ, edge = (2,3), std = stdᵖ)
add_user!(ntwᵇⁿ, node = 4)

# return the network 
return ntwᵇⁿ