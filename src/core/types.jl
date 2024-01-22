################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# abstract types
abstract type AbstractUGF end

abstract type AbstractInfo end

abstract type AbstractDistribution{X,Y,Z} end

abstract type AbstractSTD{T} <: Graphs.AbstractGraph{T} end

abstract type AbstractNetwork{T} <: Graphs.AbstractGraph{T} end

abstract type AbstractStochasticProcess end
abstract type AbstractMarkovProcess <: AbstractStochasticProcess end
abstract type AbstractSemiMarkovProcess <: AbstractStochasticProcess end

# union types
UIE = Union{Int,Graphs.AbstractEdge}
Single = Union{Bool,Number,String,Symbol,AbstractSTD,Tuple}

# dict types
LibDict = Dict{UIE,Vector{Int}}
PropDict = Dict{Symbol,Any}

# broadcastable
Broadcast.broadcastable(dst::AbstractDistribution) = Ref(dst)
Broadcast.broadcastable(std::AbstractSTD) = Ref(std)
