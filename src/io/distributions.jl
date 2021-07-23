################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

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
    Exponential

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

rate(dst::Exponential) = 1.0 / dst.θ

minimum(dst::Exponential) = zero(dst.θ)
maximum(dst::Exponential) = (Inf)unit(dst.θ)

xv(dst::Exponential, z::Real) = z * dst.θ
quantile(dst::Exponential, p::Real)  = -xv(dst, log1p(-p))
cquantile(dst::Exponential, p::Real) = -xv(dst, log(p))
sojourn(dst::Exponential,dφ::Number,tol::Real) = 
    0.0unit(dφ):dφ:cquantile(dst,tol)

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

minimum(dst::Weibull) = zero(dst.θ)
maximum(dst::Weibull) = (Inf)unit(dst.θ)

xv(dst::Weibull, z::Real) = dst.θ * z ^ (1 / dst.α)
quantile(dst::Weibull, p::Real)  = xv(dst, -log1p(-p))
cquantile(dst::Weibull, p::Real) = xv(dst, -log(p))
sojourn(dst::Weibull,dφ::Number,tol::Real) = 0.0unit(dφ):dφ:cquantile(dst,tol) 

function pdf(dst::Weibull, x::Number)
    θ, α, ω = params(dst)
    dimension(θ)==dimension(x) || return false
    if x >= (0)unit(θ)
        y = (uconvert(unit(θ),x) + 1e-10unit(θ)) / θ
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

## raised cosine
# struct
"""
    Cosine

The [**raised cosine distribution**](https://en.wikipedia.org/wiki/Raised_cosine_distribution)
with mean `μ` [`Number`], maximal deviation `s` [`Number`] and optional weight
`ω` [`Real`] has a probability density function
```math
f(x, μ, σ, ω) = \\begin{cases}
                    0                                                                  &\\text{if:}~x < μ-s,        \\\\
                    \\frac{ω}{2s} \\big(1 + π \\cos\\big(\\frac{x - μ}{s}\\big)\big)   &\\text{if:}~μ-s ≤ x ≤ μ+s,  \\\\
                    0                                                                  &\\text{if:}~x > μ+s.
                \\end{cases}
```
# Constructors
| Full              | Abbr.         | Description                                               |
| :---------------- | :------------ | :-------------------------------------------------------- |
| `Cosine(μ,s,ω)`   | `𝑪(μ,s,ω)`    | full constructor                                           |
| `Cosine(μ,s)`     | `𝑪(μ,s)`      | constructor which defaults to `Cosine(μ,s,1.0)`         |
| `Cosine(μ)`       | `𝑪(μ)`        | constructor which defaults to `Cosine(μ,1.0,1.0)`         |
| `Cosine()`        | `𝑪()`         | empty constructor which defaults to `Cosine(1.0,1.0,1.0)` |
# Examples
```julia-repl
julia> Cosine()             # default Raised Cosine distr. with μ = 1.0, s = 1.0 and ω = 1.0
julia> 𝑪(3.0u"minute")      # Raised Cosine distr. with μ = 3.0 min, s = 1.0 min and ω = 1.0
julia> 𝑪(5.0u"yr",4.0u"d")  # Raised Cosine distr. with μ = 5.0 yr, s = 4.0 d and ω = 1.0
julia> 𝑪(10.0,0.5,0.2)      # scaled Raised Cosine distr. with μ = 10.0, σ = 0.5 and ω = 0.2
```
"""
struct Cosine{N<:Number, R<:Real} <: AbstractDistribution{N,R}
    μ::N            # mean
    s::N            # maximal deviation
    ω::R            # weight
end

# constructors
Cosine() = Cosine(1.0, 1.0, 1.0)
Cosine(μ::Number) = Cosine(μ, oneunit(μ), 1.0)
Cosine(μ::Number, s::Number) = Cosine(μ, uconvert(unit(μ), s), 1.0)

# shortened constructors
𝑪() = Cosine()
𝑪(μ::Number) = Cosine(μ)
𝑪(μ::Number, s::Number) = Cosine(μ, uconvert(unit(μ),s))
𝑪(μ::Number, s::Number, ω::Real) = Cosine(μ, uconvert(unit(μ), s), ω)

# functions
mean(dst::Cosine)   = dst.μ
dev(dst::Cosine)    = dst.s
weight(dst::Cosine) = dst.ω
params(dst::Cosine) = (dst.μ, dst.s, dst.ω)

minimum(dst::Cosine) = dst.μ-dst.s
maximum(dst::Cosine) = dst.μ+dst.s

quantile(dst::Cosine, p::Real) = quantile_bisect(dst, p)
cquantile(dst::Cosine, p::Real) = quantile(dst, 1-p)
sojourn(dst::Cosine, dφ::Number, tol::Real) = 
    dst.μ-dst.s:dφ:dst.μ+dst.s+dφ

function pdf(dst::Cosine, x::Number)
    μ, s, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(x) || return false
    if μ-s <= x <= μ+s
        y = (uconvert(unit(μ),x) - μ) / (2 * s)
        ω / (2 * s) * (1 + cospi(y))
    else
        zero(Real)/oneunit(μ)
    end
end
function cdf(dst::Cosine, x::Number)
    μ, s, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(x) || return false
    if μ-s > x
        zero(Real)
    elseif μ-s <= x <= μ+s
        y = (uconvert(unit(μ),x) - μ) / s
        ω / 2 * (1 + y + sinpi(y) / π)
    elseif x > μ+s
        ω
    end
end
function ccdf(dst::Cosine, x::Number)
    μ, s, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(x) || return false
    if μ-s > x
        ω
    elseif μ-s <= x <= μ+s
        y = (μ - uconvert(unit(μ), x)) / s
        ω / 2 * (1 + y + sinpi(y) / π)
    elseif x > μ+s
        zero(Real)
    end
end

## dirac
"""
    Dirac

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

minimum(dst::Dirac) = dst.o
maximum(dst::Dirac) = dst.o

quantile(dst::Dirac, p::Real)  = dst.o
cquantile(dst::Dirac, p::Real) = dst.o

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

## uniform
"""
    Uniform

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

minimum(dst::Uniform) = dst.a
maximum(dst::Uniform) = dst.b

quantile(dst::Uniform, p::Real)  = dst.a + p * (dst.b - dst.a)
cquantile(dst::Uniform, p::Real) = dst.b + p * (dst.a - dst.b)
sojourn(dst::Uniform,dφ::Number,tol::Real) = 0.0unit(dφ):dφ:cquantile(dst,tol) 

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