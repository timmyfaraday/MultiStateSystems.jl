################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models, often found in           #
# reliability engineering.                                                     #
# See https://github.com/timmyfaraday/MultiStateSystems.jl                     #
################################################################################
# Authors: Tom Van Acker, Glenn Emmers                                         #
################################################################################
# Changelog:                                                                   #
# v0.3.0 - init                                                                #
################################################################################

# types ########################################################################
## abstract types
mutable struct SemiMarkovProcess <: AbstractSemiMarkovProcess end

# constants ####################################################################
## properties
const semi_markov_process_props = [:renewal, :dynamic]

# functions ####################################################################
## stochastic process
""
function get_A(std::AbstractSTD, t::StepRangeLen, tol::Real)
    # A indicates that state j can be reached if the process was initially in 
    # state i and remains there until time t. The unit of A is [unit(1/t)]
    dt = step(t)
    Ns = ns(std)
    Nt = length(t)

    A = zeros(typeof(1/dt), Ns*Nt)

    for tr in transitions(std)
        dst     = get_prop(std, tr, :distr)
        fr, to  = Graphs.src(tr), Graphs.dst(tr)
        init    = get_prop(std, fr, :init)

        init > 0.0 ? _fill_A!(to, init, Ns, tol, dst, t, A) : ~ ;
    end

    return A
end
""
function _fill_A!(to::Int, init::Real, Ns::Int, tol::Real, 
                  dst::AbstractDistribution, t::StepRangeLen, A::Vector)
    dt      = step(t)
    lb, ub  = quantile(dst, tol) - dt, cquantile(dst, tol) + dt
    
    for (ni,nt) in enumerate(t) if lb <= nt <= ub
        id = Ns * (ni-1) + to

        A[id] += init * pdf(dst, nt, zero(dt))
        isnan(A[id]) ? A[id] = zero(1/dt) : ~ ;
    end end
end
""
function get_U(std::AbstractSTD, t::StepRangeLen, tol::Real)
    # U indicates that state j can be reached if the process enters state i at 
    # sojourn time φ and remains there until calendar time t. The unit of U is [-]
    Ns = ns(std)
    Nt = length(t)

    I, J, V = Int[1:Ns*Nt;], Int[1:Ns*Nt;], ones(Float64,Ns*Nt)
  
    for tr in transitions(std)
        dst     = get_prop(std, tr, :distr)
        fr, to  = Graphs.src(tr), Graphs.dst(tr)

        _fill_U!(fr, to, Ns, tol, dst, t, I, J, V)
    end

    return _SA.sparse(I, J, V)
end
""
function _fill_U!(fr::Int, to::Int, Ns::Int, tol::Real, dst::AbstractDistribution,
                  t::StepRangeLen, I::Vector, J::Vector, V::Vector)
    dt      = step(t)
    lb, ub  = quantile(dst, tol) - dt, cquantile(dst, tol) + dt

    for (ni,nt) in enumerate(t)
        φ = nt:-dt:t[1]
        τ = reverse(φ)

        for (nj,nφ) in enumerate(φ) if lb <= nφ <= ub
            push!(I, Ns * (ni-1) + to)
            push!(J, Ns * (nj-1) + fr)
            push!(V, -dt * weights(ni,nj) * (pdf(dst, nφ, τ[nj]) |> unit(1/dt)))
            isnan(V[end]) ? V[end] = 0.0 : ~ ;
        end end
    end
end
""
function set_p!(std::AbstractSTD, t::StepRangeLen, H::Vector, tol::Real)
    # p indicates the probability of being in a state i at calendar time t. The 
    # unit of p is [-].
    Ns = ns(std)
    Nt = length(t)

    for st in states(std)
        h       = map(x -> H[Ns * (x-1) + st], 1:Nt)
        init    = get_prop(std, st, :init)
        DST     = [get_prop(std, Graphs.Edge(st,nx), :distr)
                        for nx in Graphs.outneighbors(std.graph, st)]

        _set_p!(std, st, init, tol, DST, t, h)
    end
end
""
function _set_p!(std::AbstractSTD, st::Int, init::Real, tol::Real, DST::Vector, 
                 t::StepRangeLen, h::Vector)
    # p = init ⋅ ccdf(st) + ∫ h dt - ∑ᵢ ∫ h ⋅ cdf(dstᵢ) dt
    # p = init - ∑ᵢ init * cdf(dstᵢ) + ∫ h dt - ∑ᵢ ∫ h ⋅ cdf(dstᵢ) dt 
    # p = init + ∫ h dt - ∑ᵢ init ⋅ cdf(dstᵢ) + ∫ h ⋅ cdf(dstᵢ) dt
    dt = step(t)

    p = init .+ [sum(dt * weights(ni,nj) * h[nj] for nj in 1:ni) 
                                                 for (ni,nt) in enumerate(t)]
    for dst in DST 
        _fill_p!(init, dst, t, h, p) 
    end

    set_prop!(std, st, :h, h) 
    set_prop!(std, st, :prob, p)
end
""
function _fill_p!(init::Real, dst::AbstractDistribution, t::StepRangeLen, 
                  h::Vector, p::Vector)
    dt = step(t)

    if init > 0
        for (ni,nt) in enumerate(t)
            τ = t[1]:dt:nt
            
            p[ni] -= init * cdf(dst, nt, t[1])
            p[ni] -= sum(dt * weights(ni,nj) * h[nj] * cdf(dst, nt-nτ, nτ)
                            for (nj,nτ) in enumerate(τ))
        end
    else
        for (ni,nt) in enumerate(t)
            τ = t[1]:dt:nt

            p[ni] -= sum(dt * weights(ni,nj) * h[nj] * cdf(dst, nt-nτ, nτ)
                            for (nj,nτ) in enumerate(τ))
        end
    end
end
""
function solve!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
                tsim::Number=1.0u"yr", dt::Number=4u"hr", tol::Real=1e-8)
    # get the input
    dt = dt |> unit(tsim)
    t = zero(dt):dt:tsim

    # solve the problem
    U = get_U(std, t, tol)
    A = ustrip(get_A(std, t, tol))
    H = U \ A * unit(1/dt)

    set_p!(std, t, H, tol)

    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)

    # set the solved status
    set_info!(std, :solved, true)
end
""
function state_conv(dst::MultiStateSystems.AbstractDistribution, h::Vector, t::StepRangeLen, f::Int)
    h_int = _INT.cubic_spline_interpolation(t, h)
    time = 0.0 * unit(t[1]):step(t) / f:t[end]
    h_array = [h_int(nt) for nt in time]
    quant = 0.0 * unit(t[1]):step(t) / f:cquantile(dst, 1e-10)
    ccdf_array = [ccdf(dst, nt) for nt in quant]
    
    # Perform convolution
    conv_result = _DSP.conv(h_array*step(time), [0,ccdf_array...])
    
    # Truncate the convolution result to match the length of the time array
    p = conv_result[1:length(h_array)][1:f:end]
    
    return p
end