################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

"""
    EENS(usr)

Expected Energy Not Served (EENS) [MWh]

This is the amount of electricity demand of a user, measured in MWh, that is
expected not to be met in a given year.
    - EENS(usr) gives the EENS when an user's demand equals the maximal output
      of the system for that user.
"""
EENS(usr::PropDict) =
    8760u"hr"*sum((maximum(usr[:ugf].val).-usr[:ugf].val).*usr[:ugf].prb) |> u"MWh"

"""
    GRA(usr)

Generation Ratio Availability (GRA) [-]

This is the probability of at least transferring a specific percentage of the
generation, Generation Ration (GR) [-], to a user through the network.
    - GRA(usr,GR) gives the probability of transferring at least the GR to a
      user through the network.
    - GRA(usr) gives a user's GRA for a GR ranging from 0.0 to 1.0.
    - GRO(usr,GR) gives the output towards a user for a given GR
"""
GRA(usr::PropDict) = [GRA(usr,GR) for GR in 0.0:0.01:1.0]
GRA(usr::PropDict,GR::Float64) =
    sum(usr[:ugf].prb[findfirst(x -> x .>= GRO(usr,GR),usr[:ugf].val):end])
GRO(usr::PropDict,GR::Float64) = GR*maximum(usr[:ugf].val)
