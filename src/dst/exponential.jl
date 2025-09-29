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
## abstract type
abstract type AbstractExponential{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# structs ######################################################################
## struct - θ::Number, ω::Real
""
struct ExponentialNR{X<:Number, Y<:Real} <: AbstractExponential{X,Y,Function}
    θ::X            # scale 
    ω::Y            # weight: 0.0 < ω <= 1.0
end
## struct - θ::Number, ω::Function
""
struct ExponentialNF{X<:Number, Z<:Function} <: AbstractExponential{X,Real,Z}
    θ::X            # scale
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end

# constructors #################################################################
"""
    Exponential()
    Exponential(θ::Number)
    Exponential(θ::Number, ω::Real)
    Exponential(θ::Number, ω::Function)

Construct an exponential distribution with scale parameter `θ` and optional 
weight `ω`. The weight can be either a real number (0.0 < ω <= 1.0) or a 
function of time ω(t).

The exponential distribution is commonly used in reliability engineering to 
model constant failure rates.

# Arguments
- `θ::Number`: Scale parameter (mean time to failure)
- `ω::Real`: Weight parameter (0.0 < ω <= 1.0), defaults to 1.0
- `ω::Function`: Time-dependent weight function ω(t) where 0.0 < ω(t) <= 1.0

# Example
```julia-repl
julia> dst₁ = Exponential()                    # θ=1.0, ω=1.0
julia> dst₂ = Exponential(100.0u"hr")          # θ=100hr, ω=1.0
julia> dst₃ = Exponential(100.0u"hr", 0.8)     # θ=100hr, ω=0.8
julia> dst₄ = Exponential(100.0u"hr", t->0.9)  # θ=100hr, ω(t)=0.9
```
"""
Exponential() = ExponentialNR(1.0, 1.0)
Exponential(θ::X) where {X<:Number}= ExponentialNR(θ, 1.0)
Exponential(θ::X, ω::Y) where {X<:Number,Y<:Real}= ExponentialNR(θ, ω)
Exponential(θ::X, ω::Z) where {X<:Number,Z<:Function}= ExponentialNF(θ, ω)

# shortened constructors #######################################################
"""
    𝑬()
    𝑬(θ::Number)
    𝑬(θ::Number, ω::Real)
    𝑬(θ::Number, ω::Function)

Shortened constructors for exponential distributions. Equivalent to the 
`Exponential()` constructors but with Unicode notation for convenience.

# Example
```julia-repl
julia> dst₁ = 𝑬()                       # Same as Exponential()
julia> dst₂ = 𝑬(100.0u"hr")             # Same as Exponential(100.0u"hr")
julia> dst₃ = 𝑬(100.0u"hr", 0.8)        # Same as Exponential(100.0u"hr", 0.8)
```
"""
𝑬() = Exponential()
𝑬(θ::Number) = Exponential(θ)
𝑬(θ::Number, ω::Real) = Exponential(θ, ω)
𝑬(θ::Number, ω::Function) = Exponential(θ, ω)

# functions ####################################################################
## general
scale(dst::AbstractExponential)  = dst.θ
weight(dst::AbstractExponential) = dst.ω
params(dst::AbstractExponential) = (dst.θ, dst.ω)

rate(dst::AbstractExponential) = 1.0 / dst.θ

minimum(dst::AbstractExponential) = zero(dst.θ)
maximum(dst::AbstractExponential) = (Inf)unit(dst.θ)
support(dst::AbstractExponential) = (minimum(dst), maximum(dst))

## quantile
xv(dst::AbstractExponential, z::Real) = z * dst.θ
quantile(dst::AbstractExponential, p::Real)  = -xv(dst, log1p(-p))
cquantile(dst::AbstractExponential, p::Real) = -xv(dst, log(p))
sojourn(dst::AbstractExponential,dφ::Number,tol::Real) = 
    zero(dφ):dφ:cquantile(dst,tol)

## density
""
function pdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval_param(ω,t) * (1/θ) * exp(-y/θ)
    else
        zero(1/θ)
end end
""
function cdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval_param(ω,t) * (1 - exp(-y/θ))
    else
        zero(Number)
end end
""
function ccdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval_param(ω,t) * exp(-y/θ)
    else
        eval_param(ω,t)
end end