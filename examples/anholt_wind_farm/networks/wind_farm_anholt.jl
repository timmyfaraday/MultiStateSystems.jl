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

# load pkgs
using XLSX
using Plots
using Unitful
using LinearAlgebra
using MultiStateSystems

# Pkg const
const _MSS = MultiStateSystems

# Settings for a specific analysis
cbl = false
mss = false
wtr = false
number_of_clusters = 8

# Include the stochastic models for the wind turbines and cables
# include(joinpath(_MSS.BASE_DIR,"examples/anholt_wind_farm/elements/wind_turbine.jl"))
# include(joinpath(_MSS.BASE_DIR,"examples/anholt_wind_farm/elements/sea_cable.jl"))

# Feeder 1 - Nodes 30, 62:65, 76, 109:111
ntwᶠ¹ = Network()
add_user!(ntwᶠ¹, node = 1)
add_sources!(ntwᶠ¹, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ¹, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,5),(8,3),(9,8),(10,9)],
                             std = [stdᶜᵇˡ[(8.500u"km",500u"mm^2")],stdᶜᵇˡ[(1.678u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(1.464u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ¹, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,5),(8,3),(9,8),(10,9)]) ;

# Feeder 2 - Nodes 23:29, 42, 61
ntwᶠ² = Network()
add_user!(ntwᶠ², node = 1)
add_sources!(ntwᶠ², node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ², edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)],
                             std = [stdᶜᵇˡ[(4.500u"km",500u"mm^2")],stdᶜᵇˡ[(0.603u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",240u"mm^2")],stdᶜᵇˡ[(0.603u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",240u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.809u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ², edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)]) ;

# Feeder 3 - Nodes 22, 59, 60 102:108,
ntwᶠ³ = Network()
add_user!(ntwᶠ³, node = 1)
add_sources!(ntwᶠ³, node = 2:11, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ³, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,5),(11,10)],
                             std = [stdᶜᵇˡ[(4.000u"km",500u"mm^2")],stdᶜᵇˡ[(2.100u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(0.691u"km",500u"mm^2")],stdᶜᵇˡ[(1.678u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ³, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,5),(11,10)]) ;

# Feeder 4 - Nodes 19:21, 56:58, 99:101
ntwᶠ⁴ = Network()
add_user!(ntwᶠ⁴, node = 1)
add_sources!(ntwᶠ⁴, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ⁴, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,3),(8,7),(9,2),(10,9)],
                             std = [stdᶜᵇˡ[(2.000u"km",500u"mm^2")],stdᶜᵇˡ[(2.100u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(2.100u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.873u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.873u"km",150u"mm^2")],stdᶜᵇˡ[(0.691u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(1.272u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ⁴, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,3),(8,7),(9,2),(10,9)]) ;

# Feeder 5 - Nodes 17,18, 41, 54,55, 75, 96:98
ntwᶠ⁵ = Network()
add_user!(ntwᶠ⁵, node = 1)
add_sources!(ntwᶠ⁵, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ⁵, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,6),(9,4),(10,2)],
                             std = [stdᶜᵇˡ[(2.000u"km",500u"mm^2")],stdᶜᵇˡ[(1.464u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(1.095u"km",240u"mm^2")],stdᶜᵇˡ[(1.095u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(1.095u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.873u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.691u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ⁵, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,6),(9,4),(10,2)]) ;

# Feeder 6 - Nodes 15,16, 40, 52,53, 74, 92:95
ntwᶠ⁶ = Network()
add_user!(ntwᶠ⁶, node = 1)
add_sources!(ntwᶠ⁶, node = 2:11, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ⁶, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,6),(10,4),(11,2)],
                             std = [stdᶜᵇˡ[(2.000u"km",500u"mm^2")],stdᶜᵇˡ[(1.272u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(1.095u"km",240u"mm^2")],stdᶜᵇˡ[(1.272u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(1.464u"km",150u"mm^2")],stdᶜᵇˡ[(0.691u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.809u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.873u"km",150u"mm^2")],stdᶜᵇˡ[(0.691u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ⁶, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,6),(10,4),(11,2)]) ;

# Feeder 7 - Nodes 14, 39, 51, 73, 83,84, 88:91
ntwᶠ⁷ = Network()
add_user!(ntwᶠ⁷, node = 1)
add_sources!(ntwᶠ⁷, node = 2:11, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ⁷, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9),(11,6)],
                             std = [stdᶜᵇˡ[(3.000u"km",500u"mm^2")],stdᶜᵇˡ[(1.272u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(1.272u"km",500u"mm^2")],stdᶜᵇˡ[(1.464u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(1.095u"km",240u"mm^2")],stdᶜᵇˡ[(1.272u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.691u"km",150u"mm^2")],stdᶜᵇˡ[(0.691u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.691u"km",150u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ⁷, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9),(11,6)]) ;

# Feeder 8 - Nodes 13, 50, 72, 80:82, 85:87
ntwᶠ⁸ = Network()
add_user!(ntwᶠ⁸, node = 1)
add_sources!(ntwᶠ⁸, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ⁸, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)],
                             std = [stdᶜᵇˡ[(3.000u"km",500u"mm^2")],stdᶜᵇˡ[(2.100u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(1.464u"km",240u"mm^2")],stdᶜᵇˡ[(1.272u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(0.809u"km",240u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(1.464u"km",150u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.691u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ⁸, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)]) ;

# Feeder 9 - Nodes 12, 38, 49, 69:71, 77:79
ntwᶠ⁹ = Network()
add_user!(ntwᶠ⁹, node = 1)
add_sources!(ntwᶠ⁹, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ⁹, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)],
                             std = [stdᶜᵇˡ[(4.000u"km",500u"mm^2")],stdᶜᵇˡ[(1.464u"km",500u"mm^2")],
                                    stdᶜᵇˡ[(1.272u"km",240u"mm^2")],stdᶜᵇˡ[(1.678u"km",240u"mm^2")],
                                    stdᶜᵇˡ[(0.873u"km",240u"mm^2")],stdᶜᵇˡ[(0.873u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.952u"km",150u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")],
                                    stdᶜᵇˡ[(0.809u"km",150u"mm^2")]]) : 
      add_components!(ntwᶠ⁹, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)]) ;

# Feeder 10 - Nodes 11, 37, 45:48, 66:68
ntwᶠ¹⁰ = Network()
add_user!(ntwᶠ¹⁰, node = 1)
add_sources!(ntwᶠ¹⁰, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ¹⁰, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,6),(9,8),(10,9)],
                              std = [stdᶜᵇˡ[(4.500u"km",500u"mm^2")],stdᶜᵇˡ[(1.678u"km",500u"mm^2")],
                                     stdᶜᵇˡ[(1.464u"km",240u"mm^2")],stdᶜᵇˡ[(0.952u"km",240u"mm^2")],
                                     stdᶜᵇˡ[(0.952u"km",240u"mm^2")],stdᶜᵇˡ[(0.952u"km",150u"mm^2")],
                                     stdᶜᵇˡ[(0.873u"km",150u"mm^2")],stdᶜᵇˡ[(0.809u"km",150u"mm^2")],
                                     stdᶜᵇˡ[(0.809u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ¹⁰, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,6),(9,8),(10,9)]) ;      

# Feeder 11 - Nodes 10, 31:36, 43,44
ntwᶠ¹¹ = Network()
add_user!(ntwᶠ¹¹, node = 1)
add_sources!(ntwᶠ¹¹, node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ¹¹, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,6),(10,9)],
                              std = [stdᶜᵇˡ[(5.000u"km",500u"mm^2")],stdᶜᵇˡ[(1.678u"km",500u"mm^2")],
                                     stdᶜᵇˡ[(0.873u"km",240u"mm^2")],stdᶜᵇˡ[(0.873u"km",240u"mm^2")],
                                     stdᶜᵇˡ[(0.873u"km",240u"mm^2")],stdᶜᵇˡ[(0.873u"km",150u"mm^2")],
                                     stdᶜᵇˡ[(0.809u"km",150u"mm^2")],stdᶜᵇˡ[(0.691u"km",150u"mm^2")],
                                     stdᶜᵇˡ[(0.809u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ¹¹, edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,6),(10,9)]) ;

# Feeder 12 - Nodes
ntwᶠ¹² = Network()
add_user!(ntwᶠ¹², node = 1)
add_sources!(ntwᶠ¹², node = 2:10, ntw = (ntwʷᵗ,1))
cbl ? add_components!(ntwᶠ¹², edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)],
                              std = [stdᶜᵇˡ[(5.750u"km",500u"mm^2")],stdᶜᵇˡ[(0.603u"km",500u"mm^2")],
                                     stdᶜᵇˡ[(0.603u"km",240u"mm^2")],stdᶜᵇˡ[(0.603u"km",240u"mm^2")],
                                     stdᶜᵇˡ[(0.603u"km",240u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                     stdᶜᵇˡ[(0.603u"km",150u"mm^2")],stdᶜᵇˡ[(0.603u"km",150u"mm^2")],
                                     stdᶜᵇˡ[(0.603u"km",150u"mm^2")]]) :
      add_components!(ntwᶠ¹², edge = [(2,1),(3,2),(4,3),(5,4),(6,5),(7,6),(8,7),(9,8),(10,9)]) ;

# Feeders 1-6
ntwᶠ¹⁻⁶ = Network()
add_user!(ntwᶠ¹⁻⁶, node = 1)
add_sources!(ntwᶠ¹⁻⁶, node = 1, ntw = [(ntwᶠ¹,1),(ntwᶠ²,1),(ntwᶠ³,1),(ntwᶠ⁴,1),(ntwᶠ⁵,1),(ntwᶠ⁶,1)])

# Feeders 7-12
ntwᶠ⁷⁻¹² = Network()
add_user!(ntwᶠ⁷⁻¹², node = 1)
add_sources!(ntwᶠ⁷⁻¹², node = 1, ntw = [(ntwᶠ⁷,1),(ntwᶠ⁸,1),(ntwᶠ⁹,1),(ntwᶠ¹⁰,1),(ntwᶠ¹¹,1),(ntwᶠ¹²,1)])

# Overall network
ntw = Network()
add_user!(ntw, node = 1, ind = [:GRA,:EENS])
add_sources!(ntw, node = 1, ntw = [(ntwᶠ¹⁻⁶,1),(ntwᶠ⁷⁻¹²,1)])

# Solve the problem
@time solve!(ntw)

# Plot the probability distribution (p>=1e-5)
ugf = ntw.usr[1][:ugf]
scatter(ustrip.(ugf.val[ugf.prb.>=1e-5]),ugf.prb[ugf.prb.>=1e-5],legend=false)