################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# abstract types
mutable struct SemiMarkovProcess <: AbstractSemiMarkovProcess end

# properties
const semi_markov_process_props = [:renewal, :dynamic]

"""
    MultiStateSystems.set_A(std::MultiStateSystems.AbstractSTD, t::StepRangeLen)

A indicates that state j can be reached if the process was initially in state i 
and remains there until time t.

The unit of A is [unit(1/t)]
""" 
function set_A(std::AbstractSTD, t::StepRangeLen)
    dt = step(t)
    Nt = length(t)

    A = zeros(typeof(1/dt), ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        for tr in transitions(std)
            id      = ns(std) * (ni-1) + _LG.dst(tr)
            # NB: GLENN same clarification
            A[id]   += get_prop(std, _LG.src(tr), :init) * 
                            pdf(get_prop(std, tr, :distr), nt, zero(dt))
    end end 
    A[isnan.(A)] .= zero(1/dt)

    return A
end

"""
    MultiStateSystems.set_U(std::MultiStateSystems.AbstractSTD, t::StepRangeLen)

U indicates that state j can be reached if the process enters state i at time φ 
and remains there until time t

The unit of U is [-]
"""
function set_U(std::AbstractSTD, t::StepRangeLen)
    dt = step(t)
    Nt = length(t)

    U = zeros(Float64, ns(std) * Nt, ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        w = weights(ni)
        l = zero(dt):dt:nt
        for (nj,nl) in enumerate(l)
            Ψ = zeros(typeof(1/dt), ns(std), ns(std))
            for tr in transitions(std)
                # NB: GLENN additional clarification distr betreden op nϕ en blijven tot nt, φ =  nt - nl
                Ψ[_LG.dst(tr), _LG.src(tr)] = 
                    pdf(get_prop(std, tr, :distr), nt - nl, nl)
            end
            rT = (ns(std) * (ni - 1) + 1):(ns(std) * ni)
            rΦ = (ns(std) * (nj - 1) + 1):(ns(std) * nj)
            if ni == nj
                U[rT,rΦ] = Matrix(1.0_LA.I, ns(std), ns(std)) .- dt .* w[nj] .* Ψ
            else
                U[rT,rΦ] = -dt .* w[nj] .* Ψ
            end
    end end
    U[isnan.(U)] .= 0.0

    return U
end

# stochastic process
function solve!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
                tsim::Number=1.0u"yr", dt::Number=1.0u"d", tol::Real=1e-8) 
    # get the input
    t   = zero(dt):dt:tsim
    Nt  = length(t)

    # solve the problem
    Φ   = zeros(Nt, ns(std))
    H   = set_U(std, t) \ set_A(std, t)
    for st in states(std)
        for (ni,nt) in enumerate(t)
            w   = weights(ni)
            # TOM: φ <<< t, zero could be higher
            l   = zero(dt):dt:nt
            # NB: ccdf(t-l,φ) where φ = 0.0, GLENN, additional clarification
            Φ[ni,st] += get_prop(std, st, :init) * ccdf(std, st, nt, zero(dt))
            Φ[ni,st] += sum(dt .* w[nj] .* H[ns(std) * (nj-1) + st] .* 
                                ccdf(std, st, nt-nl, nl) 
                                for (nj,nl) in enumerate(l))
    end end

    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)
    set_prop!(std, states(std), :prob, [Φ[:,ns] for ns in states(std)])

    # set the solved status
    set_info!(std, :solved, true)
end

"""
    MultiStateSystems.weights(x::Int)

Determine integration weights based on extended Simpson's rule. 
w[1] and w[end] = 1/3, even weights = 4/3 and uneven weights = 2/3.
"""
function weights(x::Int)
    x==1 && return [0]
    x==2 && return [1/2, 1/2]
    x==3 && return [1/3, 4/3, 1/3]
    x==4 && return [3/8, 9/8, 9/8, 3/8]
    x==5 && return (2/45) .* [7, 32, 12, 32, 7]
    x==6 && return (5/288) .* [19, 75, 50, 50, 75, 19]
    x==7 && return (1/140) .* [41, 216, 27, 272, 27, 216, 41]
    x==8 && return (7/17280) .* [751, 3577, 1323, 2989, 2989, 1323, 3577, 751]
    
    weights             = 48 * ones(x)
    weights[1:4]        = [17, 59, 43, 49]
    weights[end-3:end]  = [49, 43, 59, 17]
    return (1/48) .* weights
end 

# # TOM: Alternative formulation for the state ccdf, using the trapping info prop
# cdf(std::AbstractSTD, ns::Int, φ::Number, t::Number) = 
#     sum(cdf(get_prop(std, _LG.Edge(ns,nx), :distr), φ, t) for nx in _LG.outneighbors(std.graph, ns))
# ccdf(std::AbstractSTD, ns::Int, φ::Number, t::Number)  = 
#     ifelse(get_prop(std, ns, :trapping), 1.0, 1.0 - cdf(std, ns, φ, t))

# check whether sum of ccdf could be suffient
ccdf(std::AbstractSTD, ns::Int, φ::Quantity, t::Quantity)  = 
    1.0 - mysum(std, ns, φ, t)

function mysum(std::AbstractSTD, ns::Int, φ::Quantity, t::Quantity) 
# Summation function to overcome summation error over empty collection.
    if _LG.outneighbors(std.graph, ns) != Int64[]
        F = sum(cdf(get_prop(std,_LG.Edge(ns,nx),:distr), φ, t)
            for nx in _LG.outneighbors(std.graph, ns))
    else
        F = 0
end end

# TOM: added functionality to enable \ with units.
elunit(A::Matrix{U}) where U = _UF.unit(U)
elunit(b::Vector{U}) where U = _UF.unit(U)
_LA.:\(A::Matrix{<:Number}, b::Vector{<:Number}) = 
    (ustrip.(A) \ ustrip(b)) * (elunit(b) / elunit(A))


