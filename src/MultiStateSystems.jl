################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

module MultiStateSystems

# import pkgs
import Base: maximum

import LightGraphs
import LinearAlgebra
import Multigraphs
import OrdinaryDiffEq

# using Pkgs
using Unitful
using AdditionalUnits

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

include("form/markov_process.jl")
include("form/steady_state_process.jl")

include("io/distributions.jl")
include("io/measure.jl")
include("io/network.jl")
include("io/state_transition_diagram.jl")
include("io/universal_generating_function.jl")
include("io/utils.jl")

include("prob/indices.jl")
include("prob/universal_generating_operator.jl")

# export
export  BASE_DIR
export  SteadyStateProcess
export  Cosine, 𝑪, Dirac, 𝑫, Exponential, 𝑬, Gamma, 𝑮, Uniform, 𝑼, Weibull, 𝑾
export  Network, add_source!, add_sources!, add_user!, add_users!,
        add_component!, add_components!
export  STD, add_state!, add_states!, add_transition!, add_transitions!
export  UGF
export  solve!

end
