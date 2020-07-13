################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

module MultiStateSystems

    # import Pkgs
    import Clustering
    import LightGraphs
    import LinearAlgebra
    import Multigraphs
    import OrdinaryDiffEq

    # using Pkgs
    using  Unitful

    # pkg const
    const _CL  = Clustering
    const _LA  = LinearAlgebra
    const _LG  = LightGraphs
    const _MG  = Multigraphs
    const _ODE = OrdinaryDiffEq
    const _UF  = Unitful

    # paths
    const BASE_DIR = dirname(@__DIR__)

    # sets
    const MsrSet = Set([:flow])

    # abstract types
    abstract type AbstractSTD{T} <: _LG.AbstractGraph{T} end

    # types
    const UIE = Union{Int,_LG.AbstractEdge}
    const LibDict = Dict{UIE,Vector{Int}}
    const PropDict = Dict{Symbol,Any}
    const Single = Union{Bool,Number,String,Symbol,AbstractSTD}

    include("io/cluster.jl")
    include("io/core.jl")
    include("io/network.jl")
    include("io/state_transition_diagram.jl")

    include("core/stochastic_process.jl")
    include("core/universal_generating_function.jl")

    # export functions
    export  BASE_DIR
    export  cluster_wind_power
    export  Network, add_source!, add_sources!, add_user!, add_users!,
            add_component!, add_components!
    export  STD, add_state!, add_states!, add_transition!, add_transitions!
    export  solve!

end
