################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Test
using Unitful
using AdditionalUnits
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

@testset "MultiStateSystems" begin
    include("indices.jl")
end
