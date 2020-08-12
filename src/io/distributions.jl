################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

## Exponential
"""
    Exponential(θ,ω)

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
# struct
struct Exponential{N<:Number, R<:Real} <: AbstractDistribution{N,R}
    θ::N            # scale
    ω::R            # weight: 0.0 < ω <= 1.0
end

# constructors
Exponential() = Exponential(1.0, 1.0)
Exponential(θ::Number) = Exponential(θ, 1.0)

# shortened constructors
𝑬() = Exponential()
𝑬(θ::Number) = Exponential(θ)
𝑬(θ::Number, ω::Real) = Exponential(θ, ω)

# functions
scale(dst::Exponential)  = dst.θ
weight(dst::Exponential) = dst.ω
params(dst::Exponential) = (dst.θ, dst.ω)
function pdf(dst::Exponential, x::Number)
    θ, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = uconvert(unit(θ),x)
        ω * (1/θ) * exp(-y/θ)
    else
        zero(R)unit(θ)
    end
end
function cdf(dst::Exponential, x::Number)
    θ, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = uconvert(unit(θ),x)
        ω * (1 - exp(-y/θ))
    else
        zero(R)
    end
end
function ccdf(dst::Exponential, x::Number)
    θ, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = uconvert(unit(θ),x)
        ω * exp(-y/θ)
    else
        ω
    end
end

## Weibull
"""
    Weibull(θ,α,ω)

The [*Weibull distribution*](http://en.wikipedia.org/wiki/Weibull_distribution)
with scale parameter `θ`, shape parameter `α` and optional weight `ω` has a
probability density function

```math
f(x; θ, α, ω) = \begin{cases}
                    (αω)/θ (x/θ)^{α-1} e^{-(x/θ)^{α}}   & x ≥ 0, \\
                    0                                   & x < 0.
                \end{cases}
```
"""
# struct
struct Weibull{N<:Number, R<:Real} <: AbstractDistribution{N,R}
    θ::N            # scale
    α::R            # shape
    ω::R            # weight: 0.0 < ω <= 1.0
end

# constructors
Weibull() = Weibull(1.0, 1.0, 1.0)
Weibull(θ::Number) = Weibull(θ, 1.0, 1.0)
Weibull(θ::Number, α::Real) = Weibull(θ, α, 1.0)

# shortened constructors
𝑾() = Weibull()
𝑾(θ::Number) = Weibull(θ)
𝑾(θ::Number, α::Real) = Weibull(θ, α)
𝑾(θ::Number, α::Real, ω::Real) = Weibull(θ, α, ω)

# functions
scale(dst::Weibull)  = dst.θ
shape(dst::Weibull)  = dst.α
weight(dst::Weibull) = dst.ω
params(dst::Weibull) = (dst.θ, dst.α, dst.ω)
function pdf(dst::Weibull, x::Number)
    θ, α, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = uconvert(unit(θ),x) / θ
        ω * (α / θ) * y^(α - 1) * exp(-y^α)
    else
        zero(R)unit(θ)
    end
end
function cdf(dst::Weibull, x::Number)
    θ, α, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = uconvert(unit(θ),x) / θ
        ustrip(ω * (1 - exp(-y^α)))
    else
        zero(R)
    end
end
function ccdf(dst::Weibull, x::Number)
    θ, α, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = uconvert(unit(θ),x) / θ
        ω * exp(-y^α)
    else
        ω
    end
end

## LogNormal
"""
    LogNormal(μ,σ,ω)

The [*log normal distribution*](http://en.wikipedia.org/wiki/Log-normal_distribution)
with mean `μ`, standard deviation `σ` and optional weight `ω` has a probability
density function

```math
f(x; μ, σ, ω) = \begin{cases}
                    ω/(√(2π)σx) e^{-(log(x)-log(μ))^{2}/(2σ²)}  & x ≥ 0, \\
                    0                                           & x < 0.
                \end{cases}
```
"""
# struct
struct LogNormal{N<:Number, R<:Real} <: AbstractDistribution{N,R}
    μ::N            # mean
    σ::N            # std
    ω::R            # weight: 0.0 < ω <= 1.0
end

# constructors
LogNormal() = LogNormal(1.0, 1.0, 1.0)
LogNormal(μ::Number) = LogNormal(μ, (1.0)unit(μ), 1.0)
LogNormal(μ::Number, σ::Number) = LogNormal(μ, uconvert(unit(μ),σ), 1.0)

# shortened constructors
𝑳𝑵() = LogNormal()
𝑳𝑵(μ::Number) = LogNormal(μ)
𝑳𝑵(μ::Number, σ::Number) = LogNormal(μ, uconvert(unit(μ),σ))
𝑳𝑵(μ::Number, σ::Number, ω::Real) = LogNormal(μ, uconvert(unit(μ),σ), ω)

# functions
mean(dst::LogNormal)  = dst.μ
stdev(dst::LogNormal)  = dst.σ
weight(dst::LogNormal) = dst.ω
params(dst::LogNormal) = (dst.μ, dst.σ, dst.ω)
function pdf(dst::LogNormal, x::Number)
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(x) || return false
    if x >= (0)unit(μ)
        y = uconvert(unit(μ),x)
        ω / (sqrt(2π) * σ * y) * exp(-(log(y) - log(μ))^2/(2 * σ^2))
    else
        zero(R)unit(μ)
    end
end
function cdf(dst::LogNormal, x::Number)
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(x) || return false
    if x >= (0)unit(μ)
        y = uconvert(unit(μ),x)
        ω/2 + ω/2 * erf((log(y) - log(μ))/(sqrt(2) * σ))
    else
        zero(R)
    end
end
function ccdf(dst::LogNormal, x::Number)
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(x) || return false
    if x >= (0)unit(μ)
        y = uconvert(unit(μ),x)
        ω/2 - ω/2 * erf((log(y) - log(μ))/(sqrt(2) * σ))
    else
        ω
    end
end

## Dirac
"""
    Dirac(o,ω)

The [*Dirac distribution*](https://en.wikipedia.org/wiki/Dirac_delta_function)
with offset `o` and optional weight `ω` has a probability density function

```math
f(x; o, ω) = \begin{cases}
                ω   & x = 0, \\
                0   & x ≠ 0.
             \end{cases}
```
"""
# struct
struct Dirac{N<:Number, R<:Real} <: AbstractDistribution{N,R}
    o::N            # offset
    ω::R            # weight: 0.0 < ω <= 1.0
end

# constructors
Dirac() = Dirac(0.0, 1.0)
Dirac(o::Number) = Dirac(o, 1.0)

# shortened constructors
𝑫() = Dirac()
𝑫(o::Number) = Dirac(o)
𝑫(o::Number, ω::Real) = Dirac(o, ω)

# functions
offset(dst::Dirac) = dst.o
weight(dst::Dirac) = dst.ω
params(dst::Dirac) = (dst.o, dst.ω)
function pdf(dst::Dirac, x::Number)
    o, ω = params(dst)
    dimension(o)==dimension(x) || return false
    if x == o
        ω
    else
        zero(R)
    end
end
function cdf(dst::Dirac, x::Number)
    o, ω = params(dst)
    dimension(o)==dimension(x) || return false
    if x > o
        ω
    else
        zero(R)
    end
end
function ccdf(dst::Dirac, x::Number)
    o, ω = params(dst)
    dimension(o)==dimension(x) || return false
    if x > o
        zero(R)
    else
        ω
    end
end

## Uniform
"""
    Uniform(a,b,ω)
The [*continuous uniform distribution*](http://en.wikipedia.org/wiki/Uniform_distribution_(continuous))
over an interval `[a, b]` and with optional weight `ω` has a probability density
function

```math
f(x; a, b, ω) = \begin{cases}
                    ω/(b - a)   & a ≤ x ≤ b, \\
                    0           & a > x || x > b.
                \end{cases}
```
"""
# struct
struct Uniform{N<:Number, R<:Real} <: Distribution{N,R}
    a::N            # start
    b::N            # end
    ω::R            # weight: 0.0 < ω <= 1.0
end

# constructors
Uniform() = Uniform(0.0, 1.0, 1.0)
Uniform(a::Number, b::Number) = Uniform(a, b, 1.0)

# shortened constructors
𝑼() = Uniform()
𝑼(a::Number, b::Number) = Uniform(a, b)
𝑼(a::Number, b::Number, ω::Real) = Uniform(a, b, ω)

# functions
fr(dst::Uniform)     = dst.a
to(dst::Uniform)     = dst.b
weight(dst::Uniform) = dst.ω
params(dst::Uniform) = (dst.a, dst.b, dst.ω)
function pdf(dst::Uniform, x::Number)
    a, b, ω = params(dst)
    dimension(a)==dimension(b)==dimension(x) || return false
    if a <= x <= b
        ω / (b - a)
    else
        zero(R)unit(a)
    end
end
function cdf(dst::Uniform, x::Number)
    a, b, ω = params(dst)
    dimension(a)==dimension(b)==dimension(x) || return false
    if a > x
        zero(R)
    elseif a <= x <= b
        y = uconvert(unit(a),x)
        ω + ω * (y - b) / (b - a)
    else
        ω
    end
end
function ccdf(dst::Uniform, x::Number)
    a, b, ω = params(dst)
    dimension(a)==dimension(b)==dimension(x) || return false
    if a > x
        ω
    elseif a <= x <= b
        y = uconvert(unit(μ),x)
        ω * (b - y) / (b - a)
    else
        zero(R)
    end
end
