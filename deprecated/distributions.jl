################################################################################
# Copyright, 2020, Tom Van Acker                                               #
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
# abstract type 
abstract type AbstractExponential{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# struct - θ::Number, ω::Real
struct ExponentialNR{X<:Number, Y<:Real} <: AbstractExponential{X,Y,Function}
    θ::X            # scale 
    ω::Y            # weight: 0.0 < ω <= 1.0
end
# constructors
Exponential() = ExponentialNR(1.0, 1.0)
Exponential(θ::X) where {X<:Number}= ExponentialNR(θ, 1.0)
Exponential(θ::X, ω::Y) where {X<:Number,Y<:Real}= ExponentialNR(θ, ω)

# struct - θ::Number, ω::Function
struct ExponentialNF{X<:Number, Z<:Function} <: AbstractExponential{X,Real,Z}
    θ::X            # scale
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end
# constructors
Exponential(θ::X, ω::Z) where {X<:Number,Z<:Function}= ExponentialNF(θ, ω)

# shortened constructors
𝑬() = Exponential()
𝑬(θ::Number) = Exponential(θ)
𝑬(θ::Number, ω::Real) = Exponential(θ, ω)
𝑬(θ::Number, ω::Function) = Exponential(θ, ω)

# general functions
scale(dst::AbstractExponential)  = dst.θ
weight(dst::AbstractExponential) = dst.ω
params(dst::AbstractExponential) = (dst.θ, dst.ω)

rate(dst::AbstractExponential) = dst.ω / dst.θ

minimum(dst::AbstractExponential) = zero(dst.θ)
maximum(dst::AbstractExponential) = (Inf)unit(dst.θ)

# quantile
xv(dst::AbstractExponential, z::Real) = z * dst.θ
quantile(dst::AbstractExponential, p::Real)  = -xv(dst, log1p(-p))
cquantile(dst::AbstractExponential, p::Real) = -xv(dst, log(p))
sojourn(dst::AbstractExponential,dφ::Number,tol::Real) = 
    zero(dφ):dφ:cquantile(dst,tol)

# density functions
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
# abstract type
abstract type AbstractWeibull{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# struct - θ::Number, α::Real, ω::Real
struct WeibullNRR{X<:Number, Y<:Real, Z<:Real} <: AbstractWeibull{X,Y,Z}
    θ::X            # scale
    α::Y            # shape
    ω::Z            # weight: 0.0 < ω <= 1.0
end
# constructors
Weibull() = WeibullNRR(1.0, 1.0, 1.0)
Weibull(θ::X) where {X<:Number} = WeibullNRR(θ, 1.0, 1.0)
Weibull(θ::X, α::Y) where {X<:Number, Y<:Real} = WeibullNRR(θ, α, 1.0)
Weibull(θ::X, α::Y, ω::Z) where {X<:Number, Y<:Real, Z<:Real} = WeibullNRR(θ, α, ω)

# struct - θ::Number, α::Real, ω::Function
struct WeibullNRF{X<:Number, Y<:Real, Z<:Function} <: AbstractWeibull{X,Y,Z}
    θ::X            # scale
    α::Y            # shape
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end
# constructors
Weibull(θ::X, α::Y, ω::Z) where {X<:Number, Y<:Real, Z<:Function} = 
    WeibullNRF(θ, α, ω)

# shortened constructors
𝑾() = Weibull()
𝑾(θ::Number) = Weibull(θ)
𝑾(θ::Number, α::Real) = Weibull(θ, α)
𝑾(θ::Number, α::Real, ω::Real) = Weibull(θ, α, ω)
𝑾(θ::Number, α::Real, ω::Function) = Weibull(θ, α, ω)

# general functions
scale(dst::AbstractWeibull)  = dst.θ
shape(dst::AbstractWeibull)  = dst.α
weight(dst::AbstractWeibull) = dst.ω
params(dst::AbstractWeibull) = (dst.θ, dst.α, dst.ω)

minimum(dst::AbstractWeibull) = zero(dst.θ)
maximum(dst::AbstractWeibull) = (Inf)unit(dst.θ)

# quantile functions
xv(dst::AbstractWeibull, z::Real) = dst.θ * z ^ (1 / dst.α)
quantile(dst::AbstractWeibull, p::Real)  = xv(dst, -log1p(-p))
cquantile(dst::AbstractWeibull, p::Real) = xv(dst, -log(p))
sojourn(dst::AbstractWeibull,dφ::Number,tol::Real) = 
    zero(dφ):dφ:cquantile(dst,tol) 

# density functions
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

## log-normal
# abstract type
abstract type AbstractLogNormal{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# struct - μ::Number, σ::Number, ω::Real
struct LogNormalNNR{X<:Number, Y<:Number, Z<:Real} <: AbstractLogNormal{X,Y,Z}
    μ::X            # mean of the corresponding normal distribution
    σ::Y            # shape of the corresponding normal distribution
    ω::Z            # weight: 0.0 < ω <= 1.0
end
# constructors
LogNormal() = LogNormalNNR(1.0, 1.0, 1.0)
LogNormal(μ::X) where {X<:Number} = LogNormalNNR(μ, 1.0unit(μ), 1.0)
LogNormal(μ::X, σ::Y) where {X<:Number, Y<:Number} = 
    LogNormalNNR(μ, uconvert(unit(μ),σ), 1.0)
LogNormal(μ::X, σ::Y, ω::Z) where {X<:Number, Y<:Number, Z<:Real} = 
    LogNormalNNR(μ, uconvert(unit(μ),σ), ω)

# struct - μ::Number, σ::Number, ω::Function
struct LogNormalNNF{X<:Number, Y<:Number, Z<:Function} <: AbstractLogNormal{X,Y,Z}
    μ::X            # mean of the corresponding normal distribution
    σ::Y            # shape of the corresponding normal distribution
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end
# constructors
LogNormal(μ::X, σ::Y, ω::Z) where {X<:Number, Y<:Number, Z<:Function}= 
    LogNormalNNF(μ, uconvert(unit(μ),σ), ω)

# shortened constructors
𝑳() = LogNormal()
𝑳(μ::Number) = LogNormal(μ)
𝑳(μ::Number, σ::Number) = LogNormal(μ, σ)
𝑳(μ::Number, σ::Number, ω::Real) = LogNormal(μ, σ, ω)
𝑳(μ::Number, σ::Number, ω::Function) = LogNormal(μ, σ, ω)

# general functions
weight(dst::AbstractLogNormal) = dst.ω
params(dst::AbstractLogNormal) = (dst.μ, dst.σ, dst.ω)

minimum(dst::AbstractLogNormal) = zero(dst.μ)
maximum(dst::AbstractLogNormal) = (Inf)unit(dst.μ)

# quantile functions
quantile(dst::AbstractLogNormal, p::Real)  = 
    exp(dst.μ + sqrt(2) * dst.σ * _SF.erfinv(2.0 * p - 1.0))
cquantile(dst::AbstractLogNormal, p::Real) = 
    exp(dst.μ + sqrt(2) * dst.σ * _SF.erfinv(2.0 * (1.0 - p) - 1.0))
sojourn(dst::AbstractLogNormal,dφ::Number,tol::Real) = 
    zero(dst.μ):uconvert(dst.μ,dφ):cquantile(dst,tol) 

# density functions
function pdf(dst::AbstractLogNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(μ)
        x = ustrip(unit(μ), σ)
        y = uconvert(unit(μ), φ)
        z = uconvert(unit(μ/μ),(log(y) - μ)^2 / (2 * σ^2))
        eval(ω,t) / (sqrt(2 * pi) * x * y) * exp(-z)
    else
        zero(1/φ)
    end
end
function cdf(dst::AbstractLogNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(μ)
        x = uconvert(unit(μ), φ)
        y = uconvert(unit(μ/μ),(log(x) - μ) / (sqrt(2) * σ))
        eval(ω,t) / 2 * _SF.erfc(-y)
    else
        zero(Number)
    end
end
function ccdf(dst::AbstractLogNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(μ)
        x = uconvert(unit(μ), φ)
        y = uconvert(unit(μ/μ),(log(x) - μ) / (sqrt(2) * σ))
        eval(ω,t) * (1 - _SF.erfc(-y) / 2)
    else
        eval(ω,t)
    end
end