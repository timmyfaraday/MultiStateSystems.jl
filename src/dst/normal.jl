################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models, often found in           #
# reliability engineering.                                                     #
# See https://github.com/timmyfaraday/MultiStateSystems.jl                     #
################################################################################
# Authors: Tom Van Acker, Glenn Emmers                                         #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# types ########################################################################
## abstract type
abstract type AbstractNormal{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# structs ######################################################################
## struct - μ::Number, σ::Number, ω::Real
""
struct NormalNNR{X<:Number, Y<:Number, Z<:Real} <: AbstractNormal{X,Y,Z}
    μ::X            # mean
    σ::Y            # standard deviation
    ω::Z            # weight: 0.0 < ω <= 1.0
end
## struct - μ::Number, σ::Number, ω::Function
""
struct NormalNNF{X<:Number, Y<:Number, Z<:Function} <: AbstractNormal{X,Y,Z}
    μ::X            # mean
    σ::Y            # standard deviation
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end

# constructors #################################################################
"""
    Normal()
    Normal(μ::Number)
    Normal(μ::Number, σ::Number)
    Normal(μ::Number, σ::Number, ω::Real)
    Normal(μ::Number, σ::Number, ω::Function)

Construct a normal (Gaussian) distribution with mean `μ`, standard deviation `σ`, 
and optional weight `ω`.

The normal distribution is commonly used in reliability engineering for modeling 
measurement errors, natural variations, and as a limiting distribution in many 
statistical applications.

# Arguments
- `μ::Number`: Mean of the distribution, defaults to 0.0
- `σ::Number`: Standard deviation of the distribution, defaults to 1.0
- `ω::Real`: Weight parameter (0.0 < ω <= 1.0), defaults to 1.0
- `ω::Function`: Time-dependent weight function ω(t) where 0.0 < ω(t) <= 1.0

# Example
```julia-repl
julia> dst₁ = Normal()                           # μ=0.0, σ=1.0, ω=1.0
julia> dst₂ = Normal(100.0u"hr")                 # μ=100hr, σ=1hr, ω=1.0
julia> dst₃ = Normal(100.0u"hr", 15.0u"hr")      # μ=100hr, σ=15hr, ω=1.0
julia> dst₄ = Normal(100.0u"hr", 15.0u"hr", 0.95) # μ=100hr, σ=15hr, ω=0.95
```
"""
Normal() = NormalNNR(0.0, 1.0, 1.0)
Normal(μ::X) where {X<:Number} = NormalNNR(μ, 1.0unit(μ), 1.0)
Normal(μ::X, σ::Y) where {X<:Number, Y<:Number} = 
    NormalNNR(μ, uconvert(unit(μ),σ), 1.0)
Normal(μ::X, σ::Y, ω::Z) where {X<:Number, Y<:Number, Z<:Real} = 
    NormalNNR(μ, uconvert(unit(μ),σ), ω)
Normal(μ::X, σ::Y, ω::Z) where {X<:Number, Y<:Number, Z<:Function}= 
    NormalNNF(μ, uconvert(unit(μ),σ), ω)

# shortened constructors #######################################################
"""
    𝑵()
    𝑵(μ::Number)
    𝑵(μ::Number, σ::Number)
    𝑵(μ::Number, σ::Number, ω::Real)
    𝑵(μ::Number, σ::Number, ω::Function)

Shortened constructors for normal distributions. Equivalent to the 
`Normal()` constructors but with Unicode notation for convenience.

# Example
```julia-repl
julia> dst₁ = 𝑵()                            # Same as Normal()
julia> dst₂ = 𝑵(100.0u"hr")                  # Same as Normal(100.0u"hr")
julia> dst₃ = 𝑵(100.0u"hr", 15.0u"hr")       # Same as Normal(100.0u"hr", 15.0u"hr")
```
"""
𝑵() = Normal()
𝑵(μ::Number) = Normal(μ)
𝑵(μ::Number, σ::Number) = Normal(μ, σ)
𝑵(μ::Number, σ::Number, ω::Real) = Normal(μ, σ, ω)
𝑵(μ::Number, σ::Number, ω::Function) = Normal(μ, σ, ω)

# functions ####################################################################
## general
mean(dst::AbstractNormal) = dst.μ
std(dst::AbstractNormal) = dst.σ
var(dst::AbstractNormal) = dst.σ^2
weight(dst::AbstractNormal) = dst.ω
params(dst::AbstractNormal) = (dst.μ, dst.σ, dst.ω)

minimum(dst::AbstractNormal) = -Inf * unit(dst.μ)
maximum(dst::AbstractNormal) = Inf * unit(dst.μ)
support(dst::AbstractNormal) = (minimum(dst), maximum(dst))

## quantile
xv(dst::AbstractNormal, z::Real) = dst.μ + z * dst.σ
quantile(dst::AbstractNormal, p::Real) = 
    xv(dst, _SF.erfinv(2*p - 1) * √2)
cquantile(dst::AbstractNormal, p::Real) = 
    xv(dst, _SF.erfinv(2*(1-p) - 1) * √2)
median(dst::AbstractNormal) = dst.μ

## density
""
function pdf(dst::AbstractNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(φ)==dimension(t) || return false
    y = uconvert(unit(μ), φ)
    z = (y - μ) / σ
    eval_param(ω, t) * (1 / (σ * √(2π))) * exp(-z^2 / 2)
end

""
function cdf(dst::AbstractNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(φ)==dimension(t) || return false
    y = uconvert(unit(μ), φ)
    z = (y - μ) / σ
    eval_param(ω, t) * (1 + _SF.erf(z / √2)) / 2
end

""
function ccdf(dst::AbstractNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(φ)==dimension(t) || return false
    y = uconvert(unit(μ), φ)
    z = (y - μ) / σ
    eval_param(ω, t) * (1 - _SF.erf(z / √2)) / 2
end
