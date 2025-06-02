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
abstract type AbstractLogNormal{X,Y,Z} <: AbstractDistribution{X,Y,Z} end

# structs ######################################################################
## struct - μ::Number, σ::Number, ω::Real
""
struct LogNormalNNR{X<:Number, Y<:Number, Z<:Real} <: AbstractLogNormal{X,Y,Z}
    μ::X            # mean of the corresponding normal distribution
    σ::Y            # shape of the corresponding normal distribution
    ω::Z            # weight: 0.0 < ω <= 1.0
end
## struct - μ::Number, σ::Number, ω::Function
""
struct LogNormalNNF{X<:Number, Y<:Number, Z<:Function} <: AbstractLogNormal{X,Y,Z}
    μ::X            # mean of the corresponding normal distribution
    σ::Y            # shape of the corresponding normal distribution
    ω::Z            # weight: 0.0 < ω(t) <= 1.0
end

# constructors #################################################################
LogNormal() = LogNormalNNR(1.0, 1.0, 1.0)
LogNormal(μ::X) where {X<:Number} = LogNormalNNR(μ, 1.0unit(μ), 1.0)
LogNormal(μ::X, σ::Y) where {X<:Number, Y<:Number} = 
    LogNormalNNR(μ, uconvert(unit(μ),σ), 1.0)
LogNormal(μ::X, σ::Y, ω::Z) where {X<:Number, Y<:Number, Z<:Real} = 
    LogNormalNNR(μ, uconvert(unit(μ),σ), ω)
LogNormal(μ::X, σ::Y, ω::Z) where {X<:Number, Y<:Number, Z<:Function}= 
    LogNormalNNF(μ, uconvert(unit(μ),σ), ω)

# shortened constructors #######################################################
𝑳() = LogNormal()
𝑳(μ::Number) = LogNormal(μ)
𝑳(μ::Number, σ::Number) = LogNormal(μ, σ)
𝑳(μ::Number, σ::Number, ω::Real) = LogNormal(μ, σ, ω)
𝑳(μ::Number, σ::Number, ω::Function) = LogNormal(μ, σ, ω)

# functions ####################################################################
## general
mean(dst::AbstractLogNormal)   = dst.μ
std(dst::AbstractLogNormal)    = dst.σ
weight(dst::AbstractLogNormal) = dst.ω
params(dst::AbstractLogNormal) = (dst.μ, dst.σ, dst.ω)

minimum(dst::AbstractLogNormal) = zero(dst.μ)
maximum(dst::AbstractLogNormal) = (Inf)unit(dst.μ)
support(dst::AbstractLogNormal) = (minimum(dst), maximum(dst))

## quantile
quantile(dst::AbstractLogNormal, p::Real)  = 
    exp(dst.μ + sqrt(2) * dst.σ * _SF.erfinv(2.0 * p - 1.0))
cquantile(dst::AbstractLogNormal, p::Real) = 
    exp(dst.μ + sqrt(2) * dst.σ * _SF.erfinv(2.0 * (1.0 - p) - 1.0))
sojourn(dst::AbstractLogNormal,dφ::Number,tol::Real) = 
    zero(dst.μ):uconvert(dst.μ,dφ):cquantile(dst,tol) 

## density
""
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
end end
""
function cdf(dst::AbstractLogNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(μ)
        x = uconvert(unit(μ), φ)
        y = uconvert(unit(μ/μ),(log(x) - μ) / (sqrt(2) * σ))
        eval(ω,t) / 2 * _SF.erfc(-y)
    else
        zero(Number)
end end
""
function ccdf(dst::AbstractLogNormal, φ::Number, t::Number=zero(φ))
    μ, σ, ω = params(dst)
    dimension(μ)==dimension(σ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(μ)
        x = uconvert(unit(μ), φ)
        y = uconvert(unit(μ/μ),(log(x) - μ) / (sqrt(2) * σ))
        eval(ω,t) * (1 - _SF.erfc(-y) / 2)
    else
        eval(ω,t)
end end