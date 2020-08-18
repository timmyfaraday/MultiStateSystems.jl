################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

## Exponential
# struct
"""
    Exponential

The [**exponential distribution**](http://en.wikipedia.org/wiki/Exponential_distribution)
with scale parameter `θ` [`Number`] and an optional weight `ω` [`Real`] has a
probability density function
```math
f(x, θ, ω) = \\begin{cases}
                \\frac{ω}{θ} \\cdot e^{-\\frac{x}{θ}}   &\\text{if:}~x ≥ 0, \\\\
                0                                       &\\text{if:}~x < 0.
             \\end{cases}
```
# Constructors
| Full              | Abbr.         | Description                                               |
| :---------------- | :------------ | :-------------------------------------------------------- |
| `Exponential(θ,ω)`| `𝑬(θ,ω)`     | full constructor                                           |
| `Exponential(θ)`  | `𝑬(θ)`       | constructor which defaults to `Exponential(θ,1.0)`         |
| `Exponential()`   | `𝑬()`        | empty constructor which defaults to `Exponential(1.0,1.0)` |
# Examples
```julia-repl
julia> Exponential()        # default exp. distr. with θ = 1.0 and ω = 1.0
julia> 𝑬(3.0u"yr")          # exp. distr. with θ = 3.0 yr and ω = 1.0
julia> 𝑬(1.0u"hr",0.2)      # scaled exp. distr. with θ = 1.0 hr and ω = 0.2
```
"""
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
# struct
"""
    Weibull

The [**Weibull distribution**](http://en.wikipedia.org/wiki/Weibull_distribution)
with scale parameter `θ` [`Number`], shape parameter `α` [`Real`] and optional
weight `ω` [`Real`] has a probability density function
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
| `Weibull(θ,α)`    | `𝑾(θ,α)`     | constructor which defaults to `Weibull(θ,α,1.0)`         |
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
# struct
"""
    LogNormal

The [**log normal distribution**](http://en.wikipedia.org/wiki/Log-normal_distribution)
with mean `μ` [`Number`], standard deviation `σ` [`Number`] and optional weight
`ω` [`Real`] has a probability density function
```math
f(x, μ, σ, ω) = \\begin{cases}
                    \\frac{ω}{\\sqrt{2π}σx} e^{-\\frac{(\\log(x)-\\log(μ))^{2}}{2σ²}}  &\\text{if:}~x ≥ 0, \\\\
                    0                                                                  &\\text{if:}~x < 0.
                \\end{cases}
```
# Constructors
| Full              | Abbr.         | Description                                               |
| :---------------- | :------------ | :-------------------------------------------------------- |
| `LogNormal(μ,σ,ω)`| `𝑳𝑵(μ,σ,ω)`  | full constructor                                           |
| `LogNormal(μ,σ)`  | `𝑳𝑵(μ,σ)`    | constructor which defaults to `LogNormal(μ,σ,1.0)`         |
| `LogNormal(μ)`    | `𝑳𝑵(μ)`      | constructor which defaults to `LogNormal(μ,1.0,1.0)`         |
| `LogNormal()`     | `𝑳𝑵()`       | empty constructor which defaults to `LogNormal(1.0,1.0,1.0)` |
# Examples
```julia-repl
julia> LogNormal()           # default LogNormal distr. with μ = 1.0, σ = 1.0 and ω = 1.0
julia> 𝑳𝑵(3.0u"minute")     # LogNormal distr. with μ = 3.0 min, σ = 1.0 and ω = 1.0
julia> 𝑳𝑵(5.0u"yr",4.0)     # LogNormal distr. with μ = 5.0 yr, σ = 4.0 and ω = 1.0
julia> 𝑳𝑵(10.0,0.5,0.2)     # scaled LogNormal distr. with μ = 10.0, σ = 0.5 and ω = 0.2
```
"""
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
# struct
"""
    Dirac

The [**Dirac distribution**](https://en.wikipedia.org/wiki/Dirac_delta_function)
with offset `o` [`Number`] and optional weight `ω` [`Real`] has a probability
density function

```math
f(x, o, ω) = \\begin{cases}
                ω   &\\text{if:}~x = 0, \\\\
                0   &\\text{if:}~x ≠ 0.
              \\end{cases}
```
# Constructors
| Full          | Abbr.         | Description                                           |
| :------------ | :------------ | :---------------------------------------------------- |
| `Dirac(o,ω)`  | `𝑫(o,ω)`     | full constructor                                      |
| `Dirac(o)`    |`𝑫(o)`        | constructor which defaults to `Dirac(o,1.0)`          |
| `Dirac()`     |`𝑫()`         | empty constructor which defaults to `Dirac(0.0,1.0)`  |
# Examples
```julia-repl
julia> Dirac()              # default dirac distr. with o = 0.0 and ω = 1.0
julia> 𝑫(20.0u"hr")         # dirac distr. with o = 20.0 hr and ω = 1.0
julia> Dirac(1.0,0.5)       # scaled dirac distr. with o = 1.0 and ω = 0.5
```
"""
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
# struct
"""
    Uniform
The [**continuous uniform distribution**](http://en.wikipedia.org/wiki/Uniform_distribution_(continuous))
over an interval `[a, b]` [`Number`] and with optional weight `ω` [`Real`] has a
probability density function
```math
f(x, a, b, ω) = \\begin{cases}
                    \\frac{ω}{b - a}   & \\text{if:}~a ≤ x ≤ b,        \\\\
                    0                  & \\text{if:}~a > x~||~b < x.
                \\end{cases}
```
# Constructors
| Full              | Abbr.         | Description                                                   |
| :------------     | :------------ | :------------------------------------------------------------ |
| `Uniform(a,b,ω)`  | `𝑼(a,b,ω)`    | full constructor                                             |
| `Uniform(a,b)`    | `𝑼(a,b)`      | constructor which defaults to `Uniform(a,b,1.0)`             |
| `Uniform()`       | `𝑼()`         | empty constructor which defaults to `Uniform(0.0,1.0,1.0)`   |
# Examples
```julia-repl
julia> 𝑼()                     # default uniform distr. with a = 0.0, b = 1.0 and ω = 1.0
julia> 𝑼(1.0u"yr",2.0u"yr")    # uniform distr. with a = 1.0 yr, b = 1.0 yr and ω = 1.0
julia> Uniform(1.0,3.0,0.4)    # scaled uniform distr. with a = 1.0, b = 3.0 and ω = 0.4
```
"""
struct Uniform{N<:Number, R<:Real} <: AbstractDistribution{N,R}
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
