#  Copyright 2020, Tom Van Acker
################################################################################
# MultiStateSystems.jl
# A Julia package to solve multi-state system models.
# See http://github.com/timmyfaraday/MultiStateSystems.jl
################################################################################

# Stochastic Process
"""
    solve!(std::MultiStateSystems.AbstractSTD, tsim::Number; alg::Symbol=:nothing)

This function determines the state probabilities of the state-transition diagram
`std` for a given simulation horizon `tsim`.

Optionally, the prefered stochastic process may be provided through the named
argument `alg`, otherwise the appropriate stochastic process is determined using
the properties of the state-transition diagram.

# Example
```julia-repl
julia> solve!(stdᵍᵉⁿ, 1000u"hr", alg = :markov_process)
```
"""
function solve!(std::AbstractSTD, tsim::Number; alg::Symbol=:nothing)
    # TODO - update_std_info!(std)
    if is_a_markov_process(std,alg) solve_markov_process!(std, tsim) end
    set_info!(std,:solved,true)
end

## Markov Process
const markov_props = [:renewal,:markovian]
is_a_markov_process(std::AbstractSTD,alg::Symbol) =
    alg==:markov_process || prod(getfield(std.props[:info],prop) for prop in markov_props)

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
end
