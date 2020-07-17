################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
# Contributor: Gayan Abeynayake                                                #
################################################################################

"""
# Expected Energy Not Supplied
"""
EENS(usr::PropDict) =
    8760u"hr"*sum((maximum(usr[:ugf].val).-usr[:ugf].val).*usr[:ugf].prb)
"""
# Generation Ratio Availability
"""
GRA(usr::PropDict,pct::Float64) =
    sum(usr[:ugf].prb[findfirst(x -> x .>= pct*maximum(usr[:ugf].val),usr[:ugf].val):end])
GRA(usr::PropDict) =
    [GRA(usr,pct) for pct in 0.0:0.01:1.0]
