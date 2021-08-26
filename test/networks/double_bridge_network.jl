################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
Example considering a double bridge network.
"""

# load pkgs 
using Unitful
using MultiStateSystems

# include the state-transition diagrams for the pipes
stdᵖ = solvedSTD(prob = [0.1, 0.9], flow = [0.0u"m^3/hr", 1.0u"m^3/hr"])

# initialize the network corresponding with the bridge network
ntwᵈᵇⁿ = Network() 

# add the sources, components and users to the network 
add_sources!(ntwᵈᵇⁿ, node = [1 ,2])
add_components!(ntwᵈᵇⁿ, edge = [(1,3), (1,4), (2,3), (2,4), (3,5), (4,5)],
                       std  = stdᵖ)
add_bidirectional_component!(ntwᵈᵇⁿ, edge = (3,4), std = stdᵖ)
add_user!(ntwᵈᵇⁿ, node = 5)

# return the network 
return ntwᵈᵇⁿ