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

# types ########################################################################
# abstract type
abstract type AbstractWeibull{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# structs ######################################################################
## struct - θ::Number, α::Real, ω::Real
""
struct WeibullNRR{X<:Number, Y<:Real, Z<:Real} <: AbstractWeibull{X,Y,Z}
    θ::X            # scale
    α::Y            # shape
    ω::Z            # weight: 0.0 < ω <= 1.0
end
## struct - θ::Number, α::Real, ω::Function
""
struct WeibullNRF{X<:Number, Y<:Real, Z<:Function} <: AbstractWeibull{X,Y,Z}
    θ::X            # scale
    α::Y            # shape
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end

# constructors #################################################################
"""
    Weibull()
    Weibull(θ::Number)
    Weibull(θ::Number, α::Real)
    Weibull(θ::Number, α::Real, ω::Real)
    Weibull(θ::Number, α::Real, ω::Function)

Construct a Weibull distribution with scale parameter `θ`, shape parameter `α`, 
and optional weight `ω`.

The Weibull distribution is widely used in reliability engineering to model 
various failure patterns, from early failures (α < 1) to wear-out failures (α > 1).

# Arguments
- `θ::Number`: Scale parameter (characteristic life), defaults to 1.0
- `α::Real`: Shape parameter, defaults to 1.0
- `ω::Real`: Weight parameter (0.0 < ω <= 1.0), defaults to 1.0
- `ω::Function`: Time-dependent weight function ω(t) where 0.0 < ω(t) <= 1.0

# Example
```julia-repl
julia> dst₁ = Weibull()                           # θ=1.0, α=1.0, ω=1.0
julia> dst₂ = Weibull(1000.0u"hr")                # θ=1000hr, α=1.0, ω=1.0
julia> dst₃ = Weibull(1000.0u"hr", 2.0)           # θ=1000hr, α=2.0, ω=1.0
julia> dst₄ = Weibull(1000.0u"hr", 2.0, 0.95)     # θ=1000hr, α=2.0, ω=0.95
```
"""
Weibull() = WeibullNRR(1.0, 1.0, 1.0)
Weibull(θ::X) where {X<:Number} = WeibullNRR(θ, 1.0, 1.0)
Weibull(θ::X, α::Y) where {X<:Number, Y<:Real} = WeibullNRR(θ, α, 1.0)
Weibull(θ::X, α::Y, ω::Z) where {X<:Number, Y<:Real, Z<:Real} = 
    WeibullNRR(θ, α, ω)
Weibull(θ::X, α::Y, ω::Z) where {X<:Number, Y<:Real, Z<:Function} = 
    WeibullNRF(θ, α, ω)

# shortened constructors #######################################################
"""
    𝑾()
    𝑾(θ::Number)
    𝑾(θ::Number, α::Real)
    𝑾(θ::Number, α::Real, ω::Real)
    𝑾(θ::Number, α::Real, ω::Function)

Shortened constructors for Weibull distributions. Equivalent to the 
`Weibull()` constructors but with Unicode notation for convenience.

# Example
```julia-repl
julia> dst₁ = 𝑾()                            # Same as Weibull()
julia> dst₂ = 𝑾(1000.0u"hr")                 # Same as Weibull(1000.0u"hr")
julia> dst₃ = 𝑾(1000.0u"hr", 2.0)            # Same as Weibull(1000.0u"hr", 2.0)
```
"""
𝑾() = Weibull()
𝑾(θ::Number) = Weibull(θ)
𝑾(θ::Number, α::Real) = Weibull(θ, α)
𝑾(θ::Number, α::Real, ω::Real) = Weibull(θ, α, ω)
𝑾(θ::Number, α::Real, ω::Function) = Weibull(θ, α, ω)

# functions ####################################################################
## general
scale(dst::AbstractWeibull)  = dst.θ
shape(dst::AbstractWeibull)  = dst.α
weight(dst::AbstractWeibull) = dst.ω
params(dst::AbstractWeibull) = (dst.θ, dst.α, dst.ω)

minimum(dst::AbstractWeibull) = zero(dst.θ)
maximum(dst::AbstractWeibull) = (Inf)unit(dst.θ)
support(dst::AbstractWeibull) = (minimum(dst), maximum(dst))

## quantile
xv(dst::AbstractWeibull, z::Real) = dst.θ * z ^ (1 / dst.α)
quantile(dst::AbstractWeibull, p::Real)  = xv(dst, -log1p(-p))
cquantile(dst::AbstractWeibull, p::Real) = xv(dst, -log(p))
sojourn(dst::AbstractWeibull,dφ::Number,tol::Real) = 
    zero(dφ):dφ:cquantile(dst,tol) 

## density 
""
function pdf(dst::AbstractWeibull, φ::Number, t::Number=zero(φ))
    θ, α, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = (uconvert(unit(θ),φ) + eps(θ)) / θ
        eval(ω,t) * (α / θ) * y^(α - 1) * exp(-y^α)
    else
        zero(1/θ)
end end
""
function cdf(dst::AbstractWeibull, φ::Number, t::Number=zero(φ))
    θ, α, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ) / θ
        ustrip(eval(ω,t) * (1 - exp(-y^α)))
    else
        zero(Number)
end end
""
function ccdf(dst::AbstractWeibull, φ::Number, t::Number=zero(φ))
    θ, α, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ) / θ
        eval(ω,t) * exp(-y^α)
    else
        eval(ω,t)
end end