################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

module MultiStateSystems

# import Pkgs
import LightGraphs
import LinearAlgebra
import Multigraphs
import OrdinaryDiffEq

# using Pkgs
using  Unitful

# pkg const
const _LA  = LinearAlgebra
const _LG  = LightGraphs
const _MG  = Multigraphs
const _ODE = OrdinaryDiffEq
const _UF  = Unitful

# paths
const BASE_DIR = dirname(@__DIR__)

# include
include("core/types.jl")

include("io/cluster.jl")
include("io/distributions.jl")
include("io/measure.jl")
include("io/network.jl")
include("io/state_transition_diagram.jl")
include("io/universal_generating_function.jl")
include("io/utils.jl")

include("prob/indices.jl")
include("prob/stochastic_process.jl")
include("prob/universal_generating_operator.jl")

# export
export  BASE_DIR
export  cluster_wind_power
export  Exponential, 𝑬, Weibull, 𝑾, LogNormal, 𝑳𝑵, Dirac, 𝑫, Uniform, 𝑼
export  pdf, cdf, ccdf  #this is for testing, may be deleted later
export  Network, add_source!, add_sources!, add_user!, add_users!,
        add_component!, add_components!
export  STD, add_state!, add_states!, add_transition!, add_transitions!
export  UGF
export  solve!

end
