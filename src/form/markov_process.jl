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
function set_rates!(std::AbstractSTD, cls::AbstractMarkovProcess)
    for nt in transitions(std)
        if !has_prop(std, nt, :rate) && has_prop(std, nt, :distr)
            if get_prop(std, nt, :distr) isa Exponential
                set_prop!(std, nt, :rate, rate(get_prop(std, nt, :distr)))
    end end end
end
function set_parameters!(std::AbstractSTD, cls::MarkovProcess)
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
function solve!(std::AbstractSTD, cls::MarkovProcess; 
                tsim::Number=1.0u"yr", dt::Number=1.0u"d", tol::Real=1e-8)
    # set the input
    set_parameters!(std, cls)

    # get the input
    t   = zero(dt):dt:tsim
    p   = [std.graph,get_tprop(std,:rate),get_sprop(std,:info)]
    u₀  = get_sprop(std,:init)
    ts  = (zero(tsim),tsim)

    # solve the problem
    if get_info(std,:time_homogeneous)
        prob = _ODE.ODEProblem(homogeneous_markov_process, u₀, ts, p)
    else
        prob = _ODE.ODEProblem(inhomogeneous_markov_process, u₀, ts, p)
    end
    sol = _ODE.solve(prob, _ODE.Tsit5(), reltol = tol, abstol = tol)

    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)
    set_prop!(std, states(std), :prob, [sol(t,idxs=ns) for ns in states(std)])

    # set the solved status
    set_info!(std, :solved, true)
end