# Stochastic Processes

## Introduction

A number of stochastic processes are available to determine the state
probabilities of an `std`:

* Steady state process `SteadyStateProcess <: AbstractMarkovProcess`
* Markov chain `:markov_process`
* Markov process `:markov_process`
* Semi-Markov process `:semimarkov_process`
* Van Acker process `:vanacker_process`

Solving a stochastic process may be accomplished through:
```@docs
solve!(std::MultiStateSystems.AbstractSTD, tsim::Number; alg::Symbol=:nothing)
```

## Steady State Process
| Spaces      |             	| Properties  |             |
| :---------- | :-------------- | :---------- | :---------- |
| State-space | discrete        | Renewal     | ✅          |
| Time-space  | singular (t=∞)	| Markov      | ✅          |

A steady state process determines the state-state probability of the state space
associated with a time-homogeneous Markov process/chain, i.e., where t→∞.

https://www.maplesoft.com/support/help/maple/view.aspx?path=examples/SteadyStateMarkovChain

## Markov Process

| Spaces      |             	| Properties  |             |
| :---------- | :-------------- | :---------- | :---------- |
| State-space | discrete    	| Renewal     | ✅          |
| Time-space  | continuous  	| Markov      | ✅          |

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

## Van Acker Process

* `T. Van Acker, and D. Van Hertem (2018). Stochastic Process for the
   Availability Assessment of Single-Feeder Industrial Energy System Sections.
   IEEE Trans. on Rel., 67(4), 1459-1467.`

| Spaces      |                 | Properties  |             |
| :---------- | :-------------- | :---------- | :---------- |
| State-space | semi-continuous | Renewal     | ❎[^1]      |
| Time-space  | continuous      | Markov      | ❎          |

[^1]: The normal operation state can be a non-renewal state, enabling imperfect
      maintenance.

An Van Acker process is described by a random variable $X_{t,φ}$, where $t$
denotes the calendar time. The possible values of $X_{t,φ}$ are represented by
the semi-continuous state-space $𝓢$ of the state transition diagram `std` and
$φ_s$ denotes its state sojourn time.

An Van Acker process is described by a single PDDE with non-local boundary
condition for the normal operation state $n ∈ 𝓝$, from which all state
probabilies of the other states may be derived.

```math
\begin{aligned}
\frac{∂p_{n}}{∂φ_{n}} + \frac{∂p_{n}}{∂t} =& - \sum_{f ∈ 𝓕} λ_{f}(t,φ_{n}) p_{n}(t,φ_{n})                            \\
                                          ~& + \sum_{f ∈ 𝓕^{\text{min}}} (𝖿_{f} * λ_{f}p_{n})(t,φ_{n})                \\
                              p_{n}(t,0)  =&   \sum_{f ∈ 𝓕^{\text{per}}} \int_{0}^{∞} (𝖿_{f} * λ_{f}p_{n})(t,φ_{n}) 𝖽φ_{f}
\end{aligned}
```

where, $𝓕$, $𝓕^{\text{min}}$ and $𝓕^{\text{per}}$ denote the overall failure
set and the failure sets respecively involving minimal and perfect maintenance.
The parameters $λ_{f}$ and $𝖿_{f}$ respectively denote the failure rate and
restoration pdf associated with a specific failure $f ∈ 𝓕$. The restoration pdf
is the sum of the convolutions of transition pdf's along all paths involving a
specific failure $f ∈ 𝓕$ but excluding that failure's transition pdf.

!!! The normal operation state is selected by setting that state's initial value
 :init to 1.0.

The problem structure permits descretization of the solution space into cohorts
$a ∈ 𝓐: t = φ_{n} + t_{a}$, where $t_{a}$ is the time for which a cohort $a$
has a zero sojourn time $φ_{n}$; translating the PDDE into an non-homogeneous
first order ODE.

```math
\begin{aligned}
\frac{𝖽p_{a,n}}{𝖽φ_{n}} =& - \sum_{f ∈ 𝓕} λ_{a,f}(φ_{n}) p_{a,n}(φ_{n})                                         \\
                        ~& + \sum_{f ∈ 𝓕^{\text{min}}} \sum_{x ∈ 𝓧} \bar{𝖿}_{f,x} λ_{a-x,f}(φ_{n}) p_{a-x,n}(φ_{n}) \\
            p_{a,n}(0)  =&   \sum_{f ∈ 𝓕^{\text{per}}} \sum_{x ∈ 𝓧} \bar{𝖿}_{f,x} λ_{a-x,f}(φ_{f}) p_{a-x,n}(φ_{f}) 𝖽φ_{f}
\end{aligned}
```

Using the solution for the normal operation state probability, all other state
probabilities $p_{p}(t,\varphi_{p}),~p ∈ 𝓟$ may be determined.

```math
\begin{aligned}
p_{p}(t,\varphi_{p})  =& \sum_{c ∈ 𝓒_{p}: f ∈ c} (𝖿^{\text{pre}}_{c,p}*p_{f})(t-φ_{p}) R_{c,p}(φ_{p})   \\
p_{f}(t)              =& \int_{0}^{\infty} λ_{f}(t,φ_{n}) p_{n}(t,\varphi_{n})𝖽φ_{n}
\end{aligned}
```
where $𝓒_{p}$ is the set of simple cycles going through the state $p$ and the
pdf $𝗳^{\text{pre}}_{c,p}$ convolutions of the transition pdf's of a cycle
starting from the normal operation state $n$ up to state $p$ excluding the
failure transition.
