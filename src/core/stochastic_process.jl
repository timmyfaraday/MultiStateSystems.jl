#  Copyright 2020, Tom Van Acker
################################################################################
# MultiStateSystems.jl
# A Julia package to solve multi-state system models.
# See http://github.com/timmyfaraday/MultiStateSystems.jl
################################################################################

"""
# Stochastic processes
"""
function solve!(std::AbstractSTD, tsim::Number)
    # TODO - update_std_info!(std)
    is_a_markov_process(std) ? solve_markov_process!(std, tsim) : ~ ;
end

"""
## Markov Process

State-space:     discrete
Time-space:      continuous
Renewal process: TRUE
Markov property: TRUE

A Markov process is described by a random variable Xₜ, where t denotes the
calendar time. The possible values of Xₜ is represented by the discrete state-
space 𝓢 of the process.

A Markov process respects the Markov property, which means it respects
    ℙ(Xₜ ∈ 𝓢 | 𝓕ₛ) = ℙ(Xₜ ∈ 𝓢 | Xₛ), ∀ s,t ∈ 𝕀: s < t,
where 𝓕ₛ represents a filtration of a probability space (Ω,𝓕,ℙ) and 𝕀 a totally
ordered index set. A Markov process is described by Kolmogorov equations, more
specifically the Kolmogorov forward equations:
    δpᵢⱼ(s;t)/δt = ∑ₖ pᵢₖ(s;t) ⋅ Aₖⱼ(t), ∀ i,j ∈ 𝓢, s,t ∈ 𝕀: s < t,
where A(t) represents the transition matrix, syn., generator matrix.

The latter describe an initial value problem for finding the state
probabilities, given transition rates ρᵢⱼ(t) and initial values δᵢ:
    dpᵢ(t)/dt = - ∑ⱼ ρᵢⱼ(t)pᵢ(t) + ∑ⱼ ρⱼᵢ(t)pⱼ(t),  ∀ i ∈ 𝓢.
"""
const markov_props = [:renewal,:markovian]
is_a_markov_process(std::AbstractSTD) =
    prod(getfield(std.props[:info],prop) for prop in markov_props)

function homogeneous_markov_process(du, u, p, t)
    G, ρ, info = p
    for ns in 1:_LG.nv(G)
        if info[ns].trapping
            du[ns] = sum(ρ[_LG.Edge(nt,ns)]*u[nt]
                            for nt in _LG.inneighbors(G,ns)
                            if ns ≠ nt)
        else
            du[ns] = sum(ρ[_LG.Edge(nt,ns)]*u[nt]
                            for nt in _LG.inneighbors(G,ns)
                            if ns ≠ nt) -
                     sum(ρ[_LG.Edge(ns,nt)]*u[ns]
                            for nt in _LG.outneighbors(G,ns)
                            if ns ≠ nt)
    end end
end
function inhomogeneous_markov_process(du, u, p, t)
    G, ρ, info = p
    for ns in 1:_LG.nv(G)
        if info[ns].trapping
            du[ns] = sum(ρ[_LG.Edge(nt,ns)](t)*u[nt]
                            for nt in _LG.inneighbors(G,ns)
                            if ns ≠ nt)
        else
            du[ns] = sum(ρ[_LG.Edge(nt,ns)](t)*u[nt]
                            for nt in _LG.inneighbors(G,ns)
                            if ns ≠ nt) -
                     sum(ρ[_LG.Edge(ns,nt)](t)*u[ns]
                            for nt in _LG.outneighbors(G,ns)
                            if ns ≠ nt)
    end end
end

function solve_markov_process!(std::STD, tsim::Number; tol::Real = 1e-8)
    p   = [std.graph,get_tprop(std,:rate),get_sprop(std,:info)]
    u₀  = get_sprop(std,:init)
    ts  = (0.0_UF.unit(tsim),tsim)

    if get_info(std,:time_homogeneous)
        prob = _ODE.ODEProblem(homogeneous_markov_process,u₀,ts,p)
    else
        prob = _ODE.ODEProblem(inhomogeneous_markov_process,u₀,ts,p)
    end
    sol = _ODE.solve(prob,_ODE.RK4(),reltol = tol,abstol = tol)

    set_prop!(std,:time,sol.t)
    set_prop!(std,states(std),:prob,[sol[ns,:] for ns in states(std)])

    set_info!(std,:solved,true)
end
