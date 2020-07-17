################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
The [Anholt Wind Farm](https://en.wikipedia.org/wiki/Anholt_Offshore_Wind_Farm)

"""

# Load Pkgs
using Unitful
using MultiStateSystems

# Include the stochastic models for the wind turbines and cables
include("../elements/wind_turbine.jl")
include("../elements/sea_cable.jl")

# Feeder 1 - Nodes 30, 62:65, 76, 109:111
ntwᶠ¹ = Network()
add_user!(ntwᶠ¹, node = 1)
add_sources!(ntwᶠ¹, node = 2:10, ntw = (ntwʷᵗ,1))
add_components!(ntwᶠ¹, edge = [(1,2),(2,3),(3,4),(4,5),(5,6),(5,7),(3,8),(8,9),(9,10)],
                       std = stdᶜᵇˡ[(10u"km",240u"mm^2")])
# Feeder 2 - Nodes 23:29, 42, 61
ntwᶠ² = Network()
add_user!(ntwᶠ², node = 1)
add_sources!(ntwᶠ², node = 2:10, ntw = (ntwʷᵗ,1))
add_components!(ntwᶠ², edge = [(1,2),(2,3),(3,4),(4,5),(5,6),(6,7),(7,8),(8,9),(9,10)],
                       std = stdᶜᵇˡ[(20u"km",240u"mm^2")])

# Feeders 1-6
ntwᶠ¹⁻⁶ = Network()
add_user!(ntwᶠ¹⁻⁶, node = 1)
add_sources!(ntwᶠ¹⁻⁶, node = 1, ntw = [(ntwᶠ¹,1),(ntwᶠ²,1),(ntwᶠ¹,1),(ntwᶠ²,1),(ntwᶠ¹,1),(ntwᶠ²,1)])

# Feeders 7-12
ntwᶠ⁷⁻¹² = Network()
add_user!(ntwᶠ⁷⁻¹², node = 1)
add_sources!(ntwᶠ⁷⁻¹², node = 1, ntw = [(ntwᶠ¹,1),(ntwᶠ²,1),(ntwᶠ¹,1),(ntwᶠ²,1),(ntwᶠ¹,1),(ntwᶠ²,1)])

# Overall network
ntw = Network()
add_user!(ntw, node = 1, ind = [:GRA,:EENS])
add_sources!(ntw, node = 1, ntw = [(ntwᶠ¹⁻⁶,1),(ntwᶠ⁷⁻¹²,1)])

# Solve the problem
solve!(ntw, type = :steady)
