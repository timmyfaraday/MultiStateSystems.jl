################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models, often found in           #
# reliability engineering.                                                     #
# See https://github.com/timmyfaraday/MultiStateSystems.jl                     #
################################################################################
# Authors: Tom Van Acker                                                       #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

module MultiStateSystems

# import pkgs
import Base: maximum

import Graphs
import Interpolations
import LinearAlgebra
import Measurements
import Multigraphs
import OrdinaryDiffEq
import QuadGK
import SparseArrays
import SpecialFunctions

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
include("core/info.jl")
include("core/msr.jl")
include("core/util.jl")

include("dst/exponential.jl")
include("dst/lognormal.jl")
include("dst/weibull.jl")

include("mss/ugf.jl")
include("mss/ugo.jl")

include("ntw/cmp.jl")
include("ntw/ntw.jl")
include("ntw/src.jl")
include("ntw/usr.jl")

include("std/state.jl")
include("std/std.jl")
include("std/trans.jl")

include("stp/markov_process.jl")
include("stp/semi_markov_process.jl")
include("stp/steady_state_process.jl")

# export
export  BASE_DIR
export  Exponential, 𝑬, LogNormal, 𝑳, Weibull, 𝑾
export  UGF
export  Network, 
        add_bidirectional_component!, add_bidirectional_components!,
        add_component!, add_components!,
        add_source!, add_sources!, 
        add_user!, add_users!
export  STD, solvedSTD, 
        add_state!, add_states!, 
        add_transition!, add_transitions!
export  SteadyStateProcess, MarkovProcess, SemiMarkovProcess
export  solve!

end
