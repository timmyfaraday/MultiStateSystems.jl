################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

## eval
eval(ω::Number,t::Number) = ω
eval(ω::Function,t::Number) = ω(t) |> u"s/s"

## quantile
quantile_bisect(dst::AbstractDistribution, p::Real) =
    quantile_bisect(dst, p, minimum(dst), maximum(dst), 1.0e-12)
function quantile_bisect(dst::AbstractDistribution, p::Real, lx::Real, rx::Real,
                         tol::Real)
    cl, cr = cdf(dst, lx), cdf(dst, rx)
    @assert cl <= p <= cr
    while rx - lx > tol
        m = (lx + rx) / 2
        c = cdf(dst, m)
        if p > c
            cl, lx = c, m
        else
            cr, rx = c, m
    end end
    return (lx + rx) / 2
end

## exponential 
"""
    MultiStateSystems.Exponential

The [*exponential distribution*](http://en.wikipedia.org/wiki/Exponential_distribution)
with scale parameter `θ` and an optional weight `ω` has a probability density
function

```math
f(x; θ, ω) = \begin{cases}
                ω/θ e^{-x/θ}    & x ≥ 0, \\
                0               & x < 0.
             \end{cases}
```
"""
# abstract type 
abstract type AbstractExponential{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# struct - θ::Number, ω::Real
struct ExponentialNR{X<:Number, Y<:Real, Z<:Function} <: AbstractExponential{X,Y,Z}
    θ::X            # scale 
    ω::Y            # weight: 0.0 < ω <= 1.0
end
# constructors
Exponential() = ExponentialNR{Number,Real,Function}(1.0, 1.0)
Exponential(θ::Number) = ExponentialNR{Number,Real,Function}(θ, 1.0)
Exponential(θ::Number, ω::Real) = ExponentialNR{Number,Real,Function}(θ, ω)

# struct - θ::Number, ω::Function
struct ExponentialNF{X<:Number, Y<:Real, Z<:Function} <: AbstractExponential{X,Y,Z}
    θ::X            # scale
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end
# constructors
Exponential(θ::Number, ω::Function) = ExponentialNF{Number,Real,Function}(θ, ω)

# shortened constructors
𝑬() = Exponential()
𝑬(θ::Number) = Exponential(θ)
𝑬(θ::Number, ω::Real) = Exponential(θ, ω)
𝑬(θ::Number, ω::Function) = Exponential(θ, ω)

# general functions
scale(dst::AbstractExponential)  = dst.θ
weight(dst::AbstractExponential) = dst.ω
params(dst::AbstractExponential) = (dst.θ, dst.ω)

rate(dst::AbstractExponential) = 1.0 / dst.θ

minimum(dst::AbstractExponential) = zero(dst.θ)
maximum(dst::AbstractExponential) = (Inf)unit(dst.θ)

xv(dst::AbstractExponential, z::Real) = z * dst.θ
quantile(dst::AbstractExponential, p::Real)  = -xv(dst, log1p(-p))
cquantile(dst::AbstractExponential, p::Real) = -xv(dst, log(p))
sojourn(dst::AbstractExponential,dφ::Number,tol::Real) = 
    zero(dφ):dφ:cquantile(dst,tol)

function pdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval(ω,t) * (1/θ) * exp(-y/θ)
    else
        zero(1/θ)
    end
end
function cdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval(ω,t) * (1 - exp(-y/θ))
    else
        zero(Number)
    end
end
function ccdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval(ω,t) * exp(-y/θ)
    else
        eval(ω,t)
    end
end

## weibull
"""
    Weibull

The [*Weibull distribution*](http://en.wikipedia.org/wiki/Weibull_distribution)
with scale parameter `θ`, shape parameter `α` and optional weight `ω` has a
probability density function

```math
f(x, θ, α, ω) = \\begin{cases}
                    \\frac{αω}{θ} \\cdot \\big(\\frac{x}{θ}\\big)^{α-1} \\cdot e^{-\\big(\\frac{x}{θ}\\big)^{α}}  &\\text{if:}~x ≥ 0, \\\\
                    0                                                                                 &\\text{if:}~x < 0.
                \\end{cases}
```
# Constructors
| Full              | Abbr.         | Description                                               |
| :---------------- | :------------ | :-------------------------------------------------------- |
| `Weibull(θ,α,ω)`  | `𝑾(θ,α,ω)`   | full constructor                                           |
| `Weibull(θ,α)`    | `𝑾(θ,α)`     | constructor which defaults to `Weibull(θ,α,1.0)`           |
| `Weibull(θ)`      | `𝑾(θ)`       | constructor which defaults to `Weibull(θ,1.0,1.0)`         |
| `Weibull()`       | `𝑾()`        | empty constructor which defaults to `Weibull(1.0,1.0,1.0)` |
# Examples
```julia-repl
julia> Weibull()            # default Weibull distr. with θ = 1.0, α = 1.0 and ω = 1.0
julia> 𝑾(3.0u"minute")     # Weibull distr. with θ = 3.0 min, α = 1.0 and ω = 1.0
julia> 𝑾(5.0u"yr",4.0)     # Weibull distr. with θ = 5.0 yr, α = 4.0 and ω = 1.0
julia> 𝑾(10.0,0.5,0.2)     # scaled Weibull distr. with θ = 10.0, α = 0.5 and ω = 0.2
```
"""
# abstract type
abstract type AbstractWeibull{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# struct - θ::Number, α::Real, ω::Real
struct WeibullNRR{X<:Number, Y<:Real, Z<:Function} <: AbstractWeibull{X,Y,Z}
    θ::X            # scale
    α::Y            # shape
    ω::Y            # weight: 0.0 < ω <= 1.0
end
# constructors
Weibull() = WeibullNRR{Number,Real,Function}(1.0, 1.0, 1.0)
Weibull(θ::Number) = WeibullNRR{Number,Real,Function}(θ, 1.0, 1.0)
Weibull(θ::Number, α::Real) = WeibullNRR{Number,Real,Function}(θ, α, 1.0)
Weibull(θ::Number, α::Real, ω::Real) = WeibullNRR{Number,Real,Function}(θ, α, ω)

# struct - θ::Number, α::Real, ω::Function
struct WeibullNRF{X<:Number, Y<:Real, Z<:Function} <: AbstractWeibull{X,Y,Z}
    θ::X            # scale
    α::Y            # shape
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end
# constructors
Weibull(θ::Number, α::Real, ω::Function) = WeibullNRF{Number,Real,Function}(θ, α, ω)

# shortened constructors
𝑾() = Weibull()
𝑾(θ::Number) = Weibull(θ)
𝑾(θ::Number, α::Real) = Weibull(θ, α)
𝑾(θ::Number, α::Real, ω::Real) = Weibull(θ, α, ω)
𝑾(θ::Number, α::Real, ω::Function) = Weibull(θ, α, ω)

# functions
scale(dst::AbstractWeibull)  = dst.θ
shape(dst::AbstractWeibull)  = dst.α
weight(dst::AbstractWeibull) = dst.ω
params(dst::AbstractWeibull) = (dst.θ, dst.α, dst.ω)

minimum(dst::AbstractWeibull) = zero(dst.θ)
maximum(dst::AbstractWeibull) = (Inf)unit(dst.θ)

xv(dst::AbstractWeibull, z::Real) = dst.θ * z ^ (1 / dst.α)
quantile(dst::AbstractWeibull, p::Real)  = xv(dst, -log1p(-p))
cquantile(dst::AbstractWeibull, p::Real) = xv(dst, -log(p))
sojourn(dst::AbstractWeibull,dφ::Number,tol::Real) = 
    zero(dφ):dφ:cquantile(dst,tol) 

function pdf(dst::AbstractWeibull, φ::Number, t::Number=zero(φ))
    θ, α, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = (uconvert(unit(θ),φ) + eps(θ)) / θ
        eval(ω,t) * (α / θ) * y^(α - 1) * exp(-y^α)
    else
        zero(1/θ)
    end
end
function cdf(dst::AbstractWeibull, φ::Number, t::Number=zero(φ))
    θ, α, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ) / θ
        ustrip(eval(ω,t) * (1 - exp(-y^α)))
    else
        zero(Number)
    end
end
function ccdf(dst::AbstractWeibull, φ::Number, t::Number=zero(φ))
    θ, α, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ) / θ
        eval(ω,t) * exp(-y^α)
    else
        eval(ω,t)
    end
end