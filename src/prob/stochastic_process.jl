################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
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
    if is_a_markov_process(std, alg) solve_markov_process!(std, tsim) end
    if is_a_vanacker_process(std, alg) solve_vanacker_process!(std, tsim) end
    set_info!(std,:solved,true)
end

## Markov Process
const markov_props = [:renewal,:markovian]
is_a_markov_process(std::AbstractSTD, alg::Symbol) =
    alg==:markov_process || prod(getfield(std.props[:info],prop) for prop in markov_props)
function set_markov_parameters!(std::AbstractSTD, tsim::Number, tol::Real)
    set_rates!(std, tol)
end
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
function solve_markov_process!(std::AbstractSTD, tsim::Number; tol::Real = 1e-8)
    set_markov_parameters!(std, tsim, tol)
    p   = [std.graph,get_tprop(std,:rate),get_sprop(std,:info)]
    u₀  = get_sprop(std,:init)
    ts  = (0.0_UF.unit(tsim),tsim)

    if get_info(std,:time_homogeneous)
        prob = _ODE.ODEProblem(homogeneous_markov_process,u₀,ts,p)
    else
        prob = _ODE.ODEProblem(inhomogeneous_markov_process,u₀,ts,p)
    end
    sol = _ODE.solve(prob,_ODE.Vern8(),reltol = tol,abstol = tol)

    set_prop!(std,:time,sol.t)
    set_prop!(std,states(std),:prob,[sol[ns,:] for ns in states(std)])
end

## Van Acker Process
const vanacker_props = []
is_a_vanacker_process(std::AbstractSTD, alg::Symbol) =
    alg == :vanacker_process || prod(getfield(std.props[:info],prop)
                                     for prop in vanacker_props)
function set_vanacker_parameters!(std::AbstractSTD, tsim::Number, tol::Real)
    # Sets
    set_cycles!(std)
    set_failures!(std)
    set_cohorts!(std, tsim, tol)
    # Parameters
    set_sojourn_time!(std, tol)
    set_convoluted_pmfs!(std, tol)
    set_failure_pmfs!(std, tol)
    set_failure_rates!(get_prop(std, :A)[1], std, tol)
end
# function solve_initial_state!(std::AbstractSTD, tsim::Number, tol::Real)
#     # Pre-allocate the necessary input
#     # Parameters
#     n       = get_prop(std, n)
#     φ       = get_prop(std, n, :φ)
#     t       = get_prop(std, :time)
#     # Booleans
#     Bᵗʰ     = get_info(std, :time_homogeneous)
#     Bᵖʳᵉᵛ   = !isempty(get_prop(std, :P))
#     Bᵐⁱⁿ    = !isempty(get_prop(std, :Fᵐⁱⁿ))
#     Bᵖᵉʳ    = !isempty(get_prop(std, :Fᵖᵉʳ)) || !isempty(get_prop(std, :Pᵖᵉʳ))
    
#     # Pre-allocation of the cohort variables
#     pa, qa  = zeros(Float64, length(φ)), zeros(Float64, length(φ))
#     Pa      = zeros(Float64, length(t))
#     A       = set_A(get_prop(std, :A)[1], std)
#     b       = set_b(get_prop(std, :A)[1], 1.0, zeros(Float64,length(φ)), std)
#     i       = OffsetArray(Float64,𝓐); i[:] = 0.0
#     nh      = zeros(Float64,length(φ), maximum(length(get_prop(std, nf, :fᵐⁱⁿ)) 
#                                                 for nf in get_prop(std, :Fᵐⁱⁿ))) 

#     # Pre-allocation of the output
#     Pn      = zeros(Float64, length(t))
#     pf      = Dict(nf => zeros(Float64, length(φ)) for nf in get_prop(std, :F))
#     Pf      = Dict(nf => zeros(Float64, length(t)) for nf in get_prop(std, :F))
#     pm      = Dict(np => 0.0 for np in get_prop(std, :P))
#     Pp      = Dict(np => zeros(Float64, length(t)) for np in get_prop(std, :P))

#     for na in get_prop(std, :A)
#         if Bᵐⁱⁿ reset_nh!(nh) end
#         if sum(b) != 0.0
#             # Solve for the cohort probability
#             pa = A \ b
#             # Account for preventive maintenance
#             if Bᵖʳᵉᵛ 
#                 adjust_p!(A, pa, qa, pm, std)
#                 interpolate_p!(na, qa, Pa, std)
#                 update_pf!(na, qa, Pa, pf, std)
#             else
#                 interpolate_p!(na, pa, Pa, std)
#                 update_pf!(na, pa, Pa, pf, std)
#             end
#             # Update the initial state probability Pn, failure probability 
#             # dictionary Pf, and preventive maintenance probability dictionary Pp
#             update_P!(na, Pa, Pn, pf, Pf, pp, Pp, tsim, std)
#             # Update the cohort variables 
#             if Bᵖᵉʳ update_i!(na, pf, pp, i, std) end
#             if Bᵐⁱⁿ update_nh!(pa, nh, std) end
#         end
#         if !Bᵗʰ
#             set_rates!(na, std, tol, trans=get_prop(std,:F))
#             update_A!(na, A, std)
#         end
#         if na != 𝓐[end] update_b!(na, i[na+1], nh[:,1], b, std)
#     end

#     set_prop!(std, n, :prob, Pn)
#     for nf in get_prop(std, :F) set_prop!(std, nf, :prob, Pf[nf]) end
#     for np in get_prop(std, :P) set_prop!(std, np, :prob, Pp[np]) end
# end
# function solve_vanacker_process!(std::AbstractSTD, tsim; tol::Real = 1e-8)
#     set_vanacker_parameters!(std, tsim, tol)
    
#     solve_initial_state!(std, tsim, tol)
# end
