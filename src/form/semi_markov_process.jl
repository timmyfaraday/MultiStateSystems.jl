################################################################################
#  Copyright 2022, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# abstract types
mutable struct SemiMarkovProcess <: AbstractSemiMarkovProcess end

# properties
const semi_markov_process_props = [:renewal, :dynamic]

"""
    MultiStateSystems.set_A(std::MultiStateSystems.AbstractSTD, t::StepRangeLen)

A indicates that state j can be reached if the process was initially in state i 
and remains there until time t.

The unit of A is [unit(1/t)]
""" 
function set_A(std::AbstractSTD, t::StepRangeLen, tol::Real)
    dt = step(t)
    Nt = length(t)

    A = zeros(typeof(1/dt), ns(std) * Nt)

    for tr in transitions(std)
        init = get_prop(std, _LG.src(tr), :init)
        if init > 0.0
            dst = get_prop(std, tr, :distr)
            idx = _LG.dst(tr)
            lb = floor(cquantile(dst, tol*1e-3) / dt) * dt
            _fill_A!(std, t, A, idx, dst, init, lb)
        end
    end
    return A
end

function _fill_A!(std, t, A, idx, dst, init, lb)
    for (ni,nt) in enumerate(t)
        if nt <= lb
            id      = ns(std) * (ni-1) + idx
            # The content of the A-vector is the probability of initially being 
            # in state (_LG.src(tr)) times the probability of transitioning out 
            # of that state entered at sojourn time zero and evaluated at 
            # time t = (nt)
            A[id]   += init * pdf(dst, nt, zero(t[1]))
            if isnan(A[id])
                A[id] = zero(1/t[1])
            end
    end end
end


"""
    MultiStateSystems.set_U(std::MultiStateSystems.AbstractSTD, t::StepRangeLen)

U indicates that state j can be reached if the process enters state i at time φ 
and remains there until time t

The unit of U is [-]
"""

function set_U(std::AbstractSTD, t::StepRangeLen, tol::Real)
    dt = step(t)
    Nt = length(t)   
    Ns = ns(std)
    # initialize three vectors I,J,V respectively row, column and value indices
    I,J,V = Int[],Int[],Float64[]

    append!(I, 1:Ns*Nt)
    append!(J, 1:Ns*Nt)
    append!(V, ones(Float64,Ns*Nt))
  
    for tr in transitions(std)
        dst = get_prop(std, tr, :distr)
        _fill_U!(tr, dt, t, Ns, dst, I, J, V, tol)
    end
    return _SA.sparse(I, J, V)
end

function _fill_U!(tr::LightGraphs.SimpleGraphs.SimpleEdge{Int64}, dt::Number, t::StepRangeLen, Ns::Int, dst::AbstractDistribution, I::Vector, J::Vector, V::Vector, tol::Real)
    lb = Int(floor(cquantile(dst, tol*1e-3) / dt))    
    for (ni,nt) in enumerate(t)
        if ni <= lb
            Φ = nt:-dt:t[1]
            # Φ represents the sojourn time
            lt = t[1]:dt:nt
            # lt is the last transition time
            NΦ = length(Φ)
            for nj in 1:NΦ
                push!(I, Ns * (ni-1) + _LG.dst(tr))
                push!(J, Ns * (nj-1) + _LG.src(tr))
                push!(V,- dt * weights(NΦ, nj) * (pdf(dst, Φ[nj], lt[nj]) |> unit(dt)^-1))
                if isnan(V[end])
                    V[end] = 0.0
                end
            end 
        else
            Φ = nt:-dt:t[1]
            # Φ represents the sojourn time
            lt = t[1]:dt:nt
            # lt is the last transition time
            NΦ = length(lt)
            for nj in NΦ:-1:NΦ-lb
                push!(I, Ns * (ni-1) + _LG.dst(tr))
                push!(J, Ns * (nj-1) + _LG.src(tr))
                push!(V,- dt * weights(NΦ, nj) * (pdf(dst, Φ[nj], lt[nj]) |> unit(dt)^-1))
                if isnan(V[end])
                    V[end] = 0.0
                end
            end 
        end
    end 
end

function solve!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
    tsim::Number=1.0u"yr", dt::Number=4u"hr", tol::Real=1e-8)
    # get the input
    dt = dt |> unit(tsim)
    t = zero(dt):dt:tsim


    # solve the problem
    U = set_U(std,t,tol)
    A = ustrip(set_A(std,t,tol))
    H=U\A*unit(1/t[1])

    set_P(std, t, H, tol)
    
    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)

    # set the solved status
    set_info!(std, :solved, true)
end


function set_P(std::AbstractSTD, t::StepRangeLen, H::Vector, tol::Real)
    Nt = length(t)

    for st in states(std)
        dst_v = [get_prop(std, _LG.Edge(st,nx),:distr) for nx in _LG.outneighbors(std.graph, st)]
        init   = get_prop(std, st, :init)
        h = map(x->H[ns(std) * (x-1) + st], 1:Nt)
        p = set_int(dst_v, t, h, init, tol)
        set_prop!(std, st, :prob, p)
    end    
end

function set_int(dst_v, t, h, init, tol)
    if init > 0.0 
        CDF = zeros(Float64, length(t))
        for dst in dst_v
            CDF = cdf!(dst, t, CDF)
        end
        p = init .* (1 .- CDF) .+ integral(dst_v, t, h, tol)
    else
        p = integral(dst_v, t, h, tol)        
    end
    return p
end


"""
    MultiStateSystems.weights(x::Int)

Determine integration weights based on Simpson's rule. 
w[1] and w[end] = 1/2, other weights = 1.
"""

function weights(N::Int, p::Int)
    if p == 1 || p == N
        return 0.5
    else
        return 1.0
    end
end

# TOM: added functionality to enable \ with units.
elunit(A::Matrix{U}) where U = _UF.unit(U)
elunit(b::Vector{U}) where U = _UF.unit(U)
_LA.:\(A::Matrix{<:Number}, b::Vector{<:Number}) = 
    (ustrip.(A) \ ustrip(b)) * (elunit(b) / elunit(A))


function integral(dst_v::Vector, t::StepRangeLen, h::Vector, tol::Real)
    d2h = abs.(diff(diff(ustrip(h))))
    id = 1
    idx = Int[1]
    reltol = maximum(ustrip(h))*tol
    # reltol = (maximum(d2h)-Base.minimum(d2h))*tol

    while true
        next_id = findfirst(x -> x > reltol, cumsum(d2h[id:end]))
        next_id == nothing ? break : ~ ;

        id += next_id
        push!(idx,id)
    end
    push!(idx,length(d2h)+2)
    τ   = t[idx]
    γ   = _INT.interpolate((τ,), h[idx], _INT.Gridded(_INT.Linear()))
    if isempty(dst_v)
        ϕ   = [_QGK.quadgk(x -> γ(x), t[1], nτ, rtol=tol)[1] for nτ in τ]
    else
        ϕ   = [_QGK.quadgk(x -> γ(x) * ccdf(dst_v, nτ-x, x), t[1], nτ, rtol=tol)[1] for nτ in τ]
    end
    p   = _INT.LinearInterpolation(τ, ϕ)(t)
    return p
end



function compress(v::Vector, steps::Int, std::AbstractSTD)
    n = length(v)
    Ns = ns(std)
    y = Int((n-Ns)/(Ns*steps))
    compressed_v = zeros((steps+1)*Ns)*unit(v[1])
    for st in 1:Ns
        for i in 1:steps+1
            compressed_v[Ns*(i-1)+st] = v[y*Ns*(i-1)+st];
        end
    end
    return compressed_v
end

ccdf(dst_v::Vector, φ::Quantity, t::Quantity) =
    1.0 - cdfsum(dst_v, φ, t)

function cdfsum(dst_v::Vector, φ::Quantity, t::Quantity)
    val = zero(Float64)
    for dst in dst_v
        val::Float64 = cdf!(dst, φ, t, val)
    end
    return val
end


function cdf!(dst::AbstractDistribution, t::StepRangeLen, C::Vector)::Vector
    C .+= cdf.(dst, t, zero(t[1]))
    return C
end


function cdf!(dst::AbstractDistribution, φ::Quantity, t::Quantity, c::Float64)::Float64
    c += cdf(dst, φ, t)
    return c
end

# function cdfsum(dst_v::Vector{T}, φ::Quantity, t::Quantity) where T<:AbstractDistribution
#     val::Float64 = 0.0
#     for dst in dst_v
#         val = cdfsum_inner(dst, φ, t, val)
#     end
#     return val
# end

# function cdfsum_inner(dst::AbstractExponential, φ::Quantity, t::Quantity, c::Float64)::Float64
#     c += cdf(dst, φ, t)
#     return c
# end

# # cdfsum_inner for the specific AbstractExponential types
# function cdfsum_inner(dst::AbstractWeibull, φ::Quantity, t::Quantity, c::Float64)::Float64
#     c += cdf(dst, φ, t)
#     return c
# end

# function cdfsum_inner(dst::AbstractLogNormal, φ::Quantity, t::Quantity, c::Float64)::Float64
#     c += cdf(dst, φ, t)
#     return c
# end

# function cdf!(dst::AbstractExponential, t::StepRangeLen, C::Vector)::Vector
#     C .+= cdf.(dst, t, zero(t[1]))
#     return C
# end

# function cdf!(dst::AbstractWeibull, t::StepRangeLen, C::Vector)::Vector
#     C .+= cdf.(dst, t, zero(t[1]))
#     return C
# end

# function cdf!(dst::AbstractLogNormal, t::StepRangeLen, C::Vector)::Vector
#     C .+= cdf.(dst, t, zero(t[1]))
#     return C
# end

function battery_system_availability(i, T, ntw_av, μ, σ)
    # Calculate the probability that the battery does not run out of charge before
    # the power sources are repaired at time t. With a lognormal distribution for 
    # the repair time with mean μ and standard deviation σ.
    p_failure = 1 - ntw_av[i]
    p_repair_time_ge_T = ccdf(LogNormal(μ, σ),T)
    return 1-p_failure*p_repair_time_ge_T
end

function battery_system_availability(i, T, ntw_av, θ)
    # Calculate the probability that the battery does not run out of charge before
    # the power sources are repaired at time t. With an exponential distribution for 
    # the repair time with mean μ and standard deviation σ.
    p_failure = 1 - ntw_av[i]
    p_repair_time_ge_T = ccdf(Exponential(θ),T)
    return 1-p_failure*p_repair_time_ge_T
end