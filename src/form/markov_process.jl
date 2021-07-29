################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# abstract types
mutable struct MarkovProcess <: AbstractMarkovProcess end

# properties
const markov_process_props = [:renewal, :markovian, :dynamic]

# parameters 
function set_parameters!(std::AbstractSTD, cls::AbstractMarkovProcess)
    set_rates!(std, cls)
end

# stochastic process 
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
function solve!(std::AbstractSTD, cls::AbstractMarkovProcess; tsim::Number=1.0u"yr", tol::Real=1e-8)
    set_parameters!(std, cls)
    p   = [std.graph,get_tprop(std,:rate),get_sprop(std,:info)]
    u₀  = get_sprop(std,:init)
    ts  = (0.0_UF.unit(tsim),tsim)

    if get_info(std,:time_homogeneous)
        prob = _ODE.ODEProblem(homogeneous_markov_process, u₀, ts, p)
    else
        prob = _ODE.ODEProblem(inhomogeneous_markov_process, u₀, ts, p)
    end
    sol = _ODE.solve(prob, _ODE.Tsit5(), reltol = tol, abstol = tol)

    set_prop!(std, :cls, cls)
    set_prop!(std, :time, sol.t)
    set_prop!(std, states(std), :prob, [sol[ns,:] for ns in states(std)])

    set_info!(std, :solved, true)
end