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
using Measurements
using AdditionalUnits
using MultiStateSystems

# pkg const
const _MSM = Measurements
const _MSS = MultiStateSystems

@testset "MultiStateSystems" begin
    # form
    include("steady_state_process.jl")
    include("markov_process.jl")

    # io
    include("indices.jl")
    include("universal_generating_function.jl")

    # prob
    include("universal_generating_operator.jl")
end
