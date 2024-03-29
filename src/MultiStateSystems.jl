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

import Graphs
import LinearAlgebra
import Multigraphs
import OrdinaryDiffEq
import Measurements
import SpecialFunctions
import QuadGK
import Interpolations
import SparseArrays

# using Pkgs
using Unitful
using AdditionalUnits

# pkg const
const _LA  = LinearAlgebra
const _MG  = Multigraphs
const _ODE = OrdinaryDiffEq
const _UF  = Unitful
const _MSM = Measurements
const _SF  = SpecialFunctions
const _QGK = QuadGK
const _INT = Interpolations
const _SA  = SparseArrays

# paths
const BASE_DIR = dirname(@__DIR__)

# include
include("core/types.jl")

include("form/markov_process.jl")
include("form/semi_markov_process.jl")
include("form/steady_state_process.jl")

include("io/distributions.jl")
include("io/indices.jl")
include("io/measure.jl")
include("io/network.jl")
include("io/state_transition_diagram.jl")
include("io/universal_generating_function.jl")
include("io/utils.jl")

include("prob/universal_generating_operator.jl")

# export
export  BASE_DIR
export  SteadyStateProcess, MarkovProcess, SemiMarkovProcess
export  Exponential, 𝑬, LogNormal, 𝑳, Weibull, 𝑾
export  Network, add_source!, add_sources!, add_user!, add_users!,
        add_component!, add_bidirectional_component!, add_components!, 
        add_bidirectional_components!
export  STD, solvedSTD, add_state!, add_states!, add_transition!, add_transitions!
export  UGF
export  solve!
export  solve_network!


end
