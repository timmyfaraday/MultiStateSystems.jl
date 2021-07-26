################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
Example considering a high voltage system with two single bus bars.


"""

# load pkgs
using Unitful
using MultiStateSystems

# include the state-transition diagrams 
stdˢ¹ = STD(prob = [0.2,0.8], power = [0.0u"MW", 1.0u"MW"])
stdˢ² = STD(prob = [0.3,0.7], power = [0.0u"MW", 1.0u"MW"])
stdᵇʳ = STD(prob = [0.1,0.9], power = [0.0u"MW", 1.0u"MW"])

# initialize the network corresponding to the high voltage system 
ntwʰᵛ = Network()

# add the users, sources and components to the network
add_users!(ntwʰᵛ, node = [5,6], eval_dep = true)
add_sources!(ntwʰᵛ, node = [1,2], name = ["source 1", "source 2"], std = [stdˢ¹, stdˢ²])
add_components!(ntwʰᵛ, edge = [(1,3),(2,4),(3,4),(3,5),(4,6)],
                       name = ["br s1","br s2", "br bb", "br u1", "br u2"],
                       std  = stdᵇʳ)

# return the network
return ntwʰᵛ