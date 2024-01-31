# Distributions

## Introduction

The package includes a collection of probabilistic distributions and related
functions. The implementation is heavily influenced by the package
[Distributions.jl](https://github.com/JuliaStats/Distributions.jl), but extends
it to allow for weighted distributions as well as the
[Unitful.jl](https://github.com/PainterQubits/Unitful.jl) parameters.

## Distributions

Each distribution has a number of constructors, both based on its full name and
abbreviation, e.g., `Weibull()` is equivalent to `𝑾()`. The abbreviation 
constructor is the first letter of the distribution in the `\bi` font.
Furthermore, additional constructors are included for limited input, where all 
remaining parameters are set to their default, e.g., `Weibull(10.0u"hr")` is 
equivalent to `Weibull(10.0u"hr",1.0,1.0)`.

All input parameters are either of the type `Number`, its subtype `Real` or a 
`Function`:
```
Number
 |- Unitful.Quantity, e.g., 10.0u"hr"
 |- Real
     |- AbstractFloat, e.g., 10.0
     |- Integer, e.g., 10
Function
```

All distributions may be scaled using a weight parameter `ω`, where
`0.0 < ω::Real ≤ 1.0` or `0.0 < ω(t)::Fuction ≤ 1.0`.

## Exponential Distribution

The [*exponential distribution*](http://en.wikipedia.org/wiki/Exponential_distribution)
with scale parameter `θ` and an optional weight `ω` has a probability density
function

```math
f(x; θ, ω) = \begin{cases}
                ω/θ e^{-x/θ}    & x ≥ 0, \\
                0               & x < 0.
             \end{cases}
```

## Weibull Distribution

The [*Weibull distribution*](http://en.wikipedia.org/wiki/Weibull_distribution)
with scale parameter `θ`, shape parameter `α` and optional weight `ω` has a
probability density function

```math
f(x, θ, α, ω) = \\begin{cases}
                    \\frac{αω}{θ} \\cdot \\big(\\frac{x}{θ}\\big)^{α-1} \\cdot e^{-\\big(\\frac{x}{θ}\\big)^{α}}  &\\text{if:}~x ≥ 0, \\\\
                    0                                                                                 &\\text{if:}~x < 0.
                \\end{cases}
```

### Constructors
| Full              | Abbr.         | Description                                               |
| :---------------- | :------------ | :-------------------------------------------------------- |
| `Weibull(θ,α,ω)`  | `𝑾(θ,α,ω)`   | full constructor                                           |
| `Weibull(θ,α)`    | `𝑾(θ,α)`     | constructor which defaults to `Weibull(θ,α,1.0)`           |
| `Weibull(θ)`      | `𝑾(θ)`       | constructor which defaults to `Weibull(θ,1.0,1.0)`         |
| `Weibull()`       | `𝑾()`        | empty constructor which defaults to `Weibull(1.0,1.0,1.0)` |

### Examples

```julia-repl
julia> Weibull()            # default Weibull distr. with θ = 1.0, α = 1.0 and ω = 1.0
julia> 𝑾(3.0u"minute")     # Weibull distr. with θ = 3.0 min, α = 1.0 and ω = 1.0
julia> 𝑾(5.0u"yr",4.0)     # Weibull distr. with θ = 5.0 yr, α = 4.0 and ω = 1.0
julia> 𝑾(10.0,0.5,0.2)     # scaled Weibull distr. with θ = 10.0, α = 0.5 and ω = 0.2
```

## LogNormal Distribution

The [*Log-normal distribution*](https://en.wikipedia.org/wiki/Log-normal_distribution)
with expected value `μ` and standard deviation `σ` of the corresponding normal
distribution and optional weight `ω` has a probability density function

```math
f(x, μ, σ, ω) = \\begin{cases}
                    \\frac{ω}{\\sqrt{2π} x σ} \\cdot \\cdot e^{-\\big(\\frac{(\\ln{x}-μ)^{2}}{2 σ^{2}}\\big)}   &\\text{if:}~x ≥ 0, \\\\
                    0                                                                                           &\\text{if:}~x < 0.
                \\end{cases}
```

Given the ln-function, all Unitful values are converted to correspond with the unit of `μ`.

### Constructors

| Full               | Abbr.         | Description                                                  |
| :----------------- | :------------ | :----------------------------------------------------------- |
| `LogNormal(μ,σ,ω)` | `𝑳(μ,σ,ω)`    | full constructor                                             |
| `LogNormal(μ,σ)`   | `𝑳(μ,σ)`      | constructor which defaults to `LogNormal(μ,σ,1.0)`           |
| `LogNormal(μ)`     | `𝑳(μ)`        | constructor which defaults to `LogNormal(μ,1.0,1.0)`         |
| `LogNormal()`      | `𝑳()`         | empty constructor which defaults to `LogNormal(1.0,1.0,1.0)` |

### Examples

```julia-repl
julia> LogNormal()          # default Log-normal distr. with μ = 1.0, σ = 1.0 and ω = 1.0
julia> 𝑳(3.0u"minute")      # Log-normal distr. with μ = 3.0 min, σ = 1.0 min and ω = 1.0
julia> 𝑳(5.0u"yr",4.0u"d")  # Log-normal distr. with μ = 5.0 yr, σ = 4.0 d and ω = 1.0
julia> 𝑳(10.0,0.5,0.2)      # scaled Log-normal distr. with μ = 10.0, σ = 0.5 and ω = 0.2
```
