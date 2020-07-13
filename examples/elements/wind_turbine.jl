################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a multi-state Markov model for a single wind turbine taken
from:

> Reliability Model of Large Wind Farms for Power System Adequacy Studies by
  A.S. Dobakhshari, and M. Fotuhi-Firuzabad (2009)

The case study in the paper considers a five-state model to represent a single
wind turbine, with the following transition matrix:

    [ 0.000 1/hr, 0.039 1/hr, 0.013 1/hr, 0.008 1/hr, 0.018 1/hr ]
    [ 0.365 1/hr, 0.000 1/hr, 0.151 1/hr, 0.045 1/hr, 0.097 1/hr ]
λ = [ 0.122 1/hr, 0.220 1/hr, 0.000 1/hr, 0.192 1/hr, 0.155 1/hr ]
    [ 0.038 1/hr, 0.093 1/hr, 0.185 1/hr, 0.000 1/hr, 0.359 1/hr ]
    [ 0.016 1/hr, 0.012 1/hr, 0.016 1/hr, 0.067 1/hr, 0.000 1/hr ]

The results do not exactly match does in the paper. However, as the underlying
data for the rates is not given, the granularity of the results can not be
improved.
"""

# Load Pkgs
using Unitful
using MultiStateSystems

# Initialize the state-transition diagrams corresponding to the output (wto) and
# reliability (wtr) of a wind turbine.
stdʷᵗᵒ = STD()
stdʷᵗʳ = STD()

# Add the states to the std's
add_states!(stdʷᵗᵒ, flow = [0.0u"MW",0.5u"MW",1.0u"MW",1.5u"MW",2.0u"MW"],
                    init = [0.0,0.0,0.0,0.0,1.0])
add_states!(stdʷᵗʳ, flow = [(Inf)u"MW",0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW",
                            0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW"],
                    init = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0])

# Add the transitions to the std's
add_transitions!(stdʷᵗᵒ, rate = [0.000u"1/hr" 0.039u"1/hr" 0.013u"1/hr" 0.008u"1/hr" 0.018u"1/hr"
                                 0.365u"1/hr" 0.000u"1/hr" 0.151u"1/hr" 0.045u"1/hr" 0.097u"1/hr"
                                 0.122u"1/hr" 0.220u"1/hr" 0.000u"1/hr" 0.192u"1/hr" 0.155u"1/hr"
                                 0.038u"1/hr" 0.093u"1/hr" 0.185u"1/hr" 0.000u"1/hr" 0.359u"1/hr"
                                 0.016u"1/hr" 0.012u"1/hr" 0.016u"1/hr" 0.067u"1/hr" 0.000u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,2),(2,1)],
                         rate = [0.0590u"1/yr",0.0132u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,3),(3,1)],
                         rate = [0.0070u"1/yr",0.1695u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,4),(4,1)],
                         rate = [0.0770u"1/yr",0.0158u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,5),(5,1)],
                         rate = [0.0420u"1/yr",0.0361u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,6),(6,1)],
                         rate = [0.0240u"1/yr",0.3704u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,7),(7,1)],
                         rate = [0.3380u"1/yr",0.0442u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,8),(8,1)],
                         rate = [0.4320u"1/yr",0.0752u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,9),(9,1)],
                         rate = [0.4370u"1/yr",0.0625u"1/hr"])
add_transitions!(stdʷᵗʳ, states = [(1,10),(10,1)],
                         rate = [0.5380u"1/yr",0.0515u"1/hr"])

# Solve the problems
solve!(stdʷᵗᵒ, 1.0u"yr")
solve!(stdʷᵗʳ, 1.0u"yr")
