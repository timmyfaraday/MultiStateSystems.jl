# Stochastic Processes

## Introduction

A number of stochastic processes are available to determine the state
probabilities of an `std`:

* Markov chain `:markov_chain`
* Markov process `:markov_process`
* Semi-Markov process `:semimarkov`
* Van Acker process `:vanacker`

Solving a stochastic process may be accomplished through:
```@docs
solve!(std::MultiStateSystems.AbstractSTD, tsim::Number; alg::Symbol=:nothing)
```

## Markov Process

| Spaces      |             | Properties  |             |
| :---------- | :---------- | :---------- | :---------: |
| State-space | discrete    | Renewal     | ✅          |
| Time-space  | continuous  | Markov      | ✅          |

A Markov process is described by a random variable $X_t$, where $t$ denotes the
calendar time. The possible values of $X_t$ are represented by the discrete
state-space 𝓢 of the state transition diagram `std`.

A Markov process respects the Markov property, which means it respects

```math
ℙ(X_t ∈ 𝓢 | 𝓕_s) = ℙ(X_t ∈ 𝓢 | X_s), ∀ s,t ∈ 𝕀: s < t,
```

where 𝓕$_s$ represents a filtration of a probability space (Ω,𝓕,ℙ) and 𝕀 a
totally ordered index set. A Markov process is described by Kolmogorov
equations, more specifically the Kolmogorov forward equations:

```math
 δp_{ij}(s;t)/δt = ∑_k p_{ik}(s;t) ⋅ A_{kj}(t), ∀ i,j ∈ 𝓢, s,t ∈ 𝕀: s < t,
```

where $A(t)$ represents the transition matrix, syn., generator matrix. The
latter may be translated into an initial value problem for finding the state
probabilities, given transition rates ρ$_{ij}$(t) and initial values δ$_{i}$:

```math
dp_i(t)/dt = - ∑_j ρ_{ij}(t)p_i(t) + ∑_j ρ_{ji}(t)p_j(t),  ∀ i ∈ 𝓢.
```
