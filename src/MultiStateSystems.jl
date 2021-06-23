################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

module MultiStateSystems

<<<<<<< HEAD
# import pkgs
import DSP
import Interpolations
=======
# import Pkgs
>>>>>>> parent of 078bbdd (Full documentation for distributions.jl)
import LightGraphs
import LinearAlgebra
import Multigraphs
import OrdinaryDiffEq
import SpecialFunctions
import Statistics

# using Pkgs
using  Unitful

# pkg const
const _DSP = DSP
const _INT = Interpolations
const _LA  = LinearAlgebra
const _LG  = LightGraphs
const _MG  = Multigraphs
const _ODE = OrdinaryDiffEq
const _SF  = SpecialFunctions
const _ST  = Statistics
const _UF  = Unitful

# paths
const BASE_DIR = dirname(@__DIR__)

# include
include("core/types.jl")

#include("io/cluster.jl")
include("io/distributions.jl")
include("io/measure.jl")
#include("io/network.jl")
include("io/state_transition_diagram.jl")
include("io/stochastic_process.jl")
#include("io/universal_generating_function.jl")
include("io/utils.jl")

include("prob/indices.jl")
include("prob/stochastic_process.jl")
#include("prob/universal_generating_operator.jl")

# export
export  BASE_DIR
export  Cosine, 𝑪, Dirac, 𝑫, Exponential, 𝑬, Gamma, 𝑮, Uniform, 𝑼, Weibull, 𝑾
export  set_vanacker_parameters! #this is for testing, may be deleted later
#export  Network, add_source!, add_sources!, add_user!, add_users!,
#        add_component!, add_components!
export  STD, add_state!, add_states!, add_transition!, add_transitions!
export  UGF
export  solve!
export  get_cycles

end
