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
Exponential() = ExponentialNR(1.0, 1.0)
Exponential(θ::X) where {X<:Number}= ExponentialNR(θ, 1.0)
Exponential(θ::X, ω::Y) where {X<:Number,Y<:Real}= ExponentialNR(θ, ω)
Exponential(θ::X, ω::Z) where {X<:Number,Z<:Function}= ExponentialNF(θ, ω)

# shortened constructors #######################################################
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
        eval(ω,t) * (1/θ) * exp(-y/θ)
    else
        zero(1/θ)
end end
""
function cdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval(ω,t) * (1 - exp(-y/θ))
    else
        zero(Number)
end end
""
function ccdf(dst::AbstractExponential, φ::Number, t::Number=zero(φ))
    θ, ω = params(dst)
    dimension(θ)==dimension(φ)==dimension(t) || return false
    if φ >= zero(θ)
        y = uconvert(unit(θ),φ)
        eval(ω,t) * exp(-y/θ)
    else
        eval(ω,t)
end end