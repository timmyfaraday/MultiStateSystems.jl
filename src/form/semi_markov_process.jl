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
function set_A(std::AbstractSTD, t::StepRangeLen)
    dt = step(t)
    Nt = length(t)

    A = zeros(typeof(1/dt), ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        for tr in transitions(std)
            id      = ns(std) * (ni-1) + _LG.dst(tr)
            # The content of the A-vector is the probability of initially being 
            # in state (_LG.src(tr)) times the probability of transitioning out 
            # of that state entered at sojourn time zero and evaluated at 
            # time t = (nt)
            A[id]   += get_prop(std, _LG.src(tr), :init) * 
                            pdf(get_prop(std, tr, :distr), nt, zero(dt))
    end end 
    A[isnan.(A)] .= zero(1/dt)

    return A
end

function set_A(std::AbstractSTD, t::StepRangeLen, H::Vector, steps::Number)
    dt = step(t)
    Nt = length(t)
    Ns = ns(std)

    A = zeros(typeof(1/dt), ns(std) * Nt)

    A[1:length(H)] = H

    for (ni,nt) in enumerate(t)
        if ni > steps + 1
            for tr in transitions(std)
                id      = ns(std) * (ni-1) + _LG.dst(tr)
                # The content of the A-vector is the probability of initially being 
                # in state (_LG.src(tr)) times the probability of transitioning out 
                # of that state entered at sojourn time zero and evaluated at 
                # time t = (nt)
                A[id]   += get_prop(std, _LG.src(tr), :init) * 
                                pdf(get_prop(std, tr, :distr), nt, zero(dt))
    end end end
    A[isnan.(A)] .= zero(1/dt)

    return A
end

"""
    MultiStateSystems.set_U(std::MultiStateSystems.AbstractSTD, t::StepRangeLen)

U indicates that state j can be reached if the process enters state i at time φ 
and remains there until time t

The unit of U is [-]
"""
# function set_U(std::AbstractSTD, t::StepRangeLen)
#     # 
#     dt = step(t)
#     Nt = length(t)

#     U = zeros(Float64, ns(std) * Nt, ns(std) * Nt)

#     for (ni,nt) in enumerate(t)
#         w = weights(ni)
#         # l = zero(dt):dt:nt
#         l = t[1]:dt:nt
#         for (nj,nl) in enumerate(l)
#             Ψ = zeros(typeof(1/dt), ns(std), ns(std))
#             for tr in transitions(std)
#                 # The results of the pdf depends on two factors:
#                     # - How long has the system been in this particular state = sojourn time (nt - nl) with nt the index of the current time t and nl the index of entrance time.
#                     # - The time at which the system entered a certain state, with nl the index of this entrance time. This determines weight attributed to this particular transition.
#                 Ψ[_LG.dst(tr), _LG.src(tr)] = 
#                     pdf(get_prop(std, tr, :distr), nt - nl, nl) |> unit(dt)^-1
#             end
#             rT = (ns(std) * (ni - 1) + 1):(ns(std) * ni)
#             rΦ = (ns(std) * (nj - 1) + 1):(ns(std) * nj)
#             if ni == nj
#                 U[rT,rΦ] = Matrix(1.0_LA.I, ns(std), ns(std)) .- dt .* w[nj] .* Ψ
#             else
#                 U[rT,rΦ] = -dt .* w[nj] .* Ψ
#             end
#     end end
#     U[isnan.(U)] .= 0.0

#     return U
# end
# function set_U(std::AbstractSTD, t::StepRangeLen)
#     # 
#     dt = step(t)
#     Nt = length(t)

#     I,J,V = Int[],Int[],Number[]

#     append!(I, 1:Ns*Nt)
#     append!(J, 1:Ns*Nt)
#     append!(V, ones(Number,Ns*Nt))

#     for (ni,nt) in enumerate(t)
#         w = weights(ni)
#         # l = zero(dt):dt:nt
#         l = t[1]:dt:nt
#         for (nj,nl) in enumerate(l)
#             Ψ = zeros(typeof(1/dt), ns(std), ns(std))
#             for tr in transitions(std)
#                 # The results of the pdf depends on two factors:
#                     # - How long has the system been in this particular state = sojourn time (nt - nl) with nt the index of the current time t and nl the index of entrance time.
#                     # - The time at which the system entered a certain state, with nl the index of this entrance time. This determines weight attributed to this particular transition.
#                 Ψ[_LG.dst(tr), _LG.src(tr)] = 
#                     pdf(get_prop(std, tr, :distr), nt - nl, nl) |> unit(dt)^-1
#             end
#             append!(I,(ns(std) * (ni - 1) + 1):(ns(std) * ni))
#             append!(J,(ns(std) * (nj - 1) + 1):(ns(std) * nj))
#             if ni == nj
#                 U[rT,rΦ] = Matrix(1.0_LA.I, ns(std), ns(std)) .- dt .* w[nj] .* Ψ
#             else
#                 U[rT,rΦ] = -dt .* w[nj] .* Ψ
#             end
#     end end
#     U[isnan.(U)] .= 0.0

#     return U
# end

# function set_U(std::AbstractSTD, t::StepRangeLen, tol::Real)
#     dt = step(t)
#     Nt = length(t)   
#     Ns = ns(std)

#     # initialize three vectors I,J,V respectively row, column and value indices
#     I,J,V = Int[],Int[],Number[]

#     append!(I, 1:Ns*Nt)
#     append!(J, 1:Ns*Nt)
#     append!(V, ones(Number,Ns*Nt))
  
#     for tr in transitions(std)
#         # dummy value not necessary, since every diagonal entry requires a one to be added to it.
#         # Sparse Arrays allow to add values to the same location twice and add those
#         # dummy_value = ifelse(_LG.src(tr) == _LG.dst(tr),zero(1/dt),oneunit(1/dt)) # dummy_value voegt nu een 1 toe wanneer dest and source gelijk zijn, maar dit moet enkel in geval ni en nj ook nog gelijk zijn. Extra controle in de tijdslus.
#         dst = get_prop(std, tr, :distr) 
#         lb = floor(cquantile(dst, tol) / dt) * dt
#         for (ni,nt) in enumerate(t)
#             Φ = nt:-dt:t[1]
#             # Φ represents the sojourn time, ranging from
#             lt = t[1]:dt:nt
#             NΦ = length(Φ)
            
#             append!(I, (Ns * (ni-1) + _LG.dst(tr)).*ones(Int,NΦ))
#             append!(J, [Ns * (nj-1) + _LG.src(tr) for nj in 1:NΦ])
#             # append!(V, .- dt .* weights(ni)[1:NΦ] .* (pdf.(dst, nt.-Φ, Φ) .|> unit(dt)^-1))
#             append!(V, .- dt .* weights(ni)[1:NΦ] .* (pdf.(dst, Φ, lt) .|> unit(dt)^-1))
#     end end
#     V[isnan.(V)] .= 0.0
#     return _SA.sparse(I, J, V)
# end


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
        # dummy value not necessary, since every diagonal entry requires a one to be added to it.
        # Sparse Arrays allow to add values to the same location twice and add those
        # dummy_value = ifelse(_LG.src(tr) == _LG.dst(tr),zero(1/dt),oneunit(1/dt)) # dummy_value voegt nu een 1 toe wanneer dest and source gelijk zijn, maar dit moet enkel in geval ni en nj ook nog gelijk zijn. Extra controle in de tijdslus.
        dst = get_prop(std, tr, :distr)
        _fill_U!(tr, dt, t, Ns, dst, I, J, V, tol)
    end
    return _SA.sparse(I, J, V)
end

function _fill_U!(tr, dt, t, Ns, dst, I, J, V, tol)
        # lb = floor(cquantile(dst, tol) / dt) * dt
    for (ni,nt) in enumerate(t)
        Φ = nt:-dt:t[1]
        # Φ represents the sojourn time, ranging from
        lt = t[1]:dt:nt
        NΦ = length(Φ)
        for nj in 1:NΦ
            push!(I, Ns * (ni-1) + _LG.dst(tr))
            push!(J, Ns * (nj-1) + _LG.src(tr))
            push!(V,- dt * weights(NΦ, nj) * (pdf(dst, Φ[nj], lt[nj]) |> unit(dt)^-1))
            if isnan(V[end])
                V[end] = 0.0
            end
        end
        # append!(V, .- dt .* weights(ni)[1:NΦ] .* (pdf.(dst, nt.-Φ, Φ) .|> unit(dt)^-1))
        # append!(V, .- dt .* weights(ni)[1:NΦ] .* (pdf.(dst, Φ, lt) .|> unit(dt)^-1))
    end 
end
   
# function set_U(std::AbstractSTD, t::StepRangeLen, steps::Number, tol::Real)
#     dt = step(t)
#     Nt = length(t)   
#     Ns = ns(std)

#     # initialize three vectors I,J,V respectively row, column and value indices
#     I,J,V = Int[],Int[],Number[]

#     append!(I, 1:Ns*Nt)
#     append!(J, 1:Ns*Nt)
#     append!(V, ones(Number,Ns*Nt))
  
#     for tr in transitions(std)
#         # dummy value not necessary, since every diagonal entry requires a one to be added to it.
#         # Sparse Arrays allow to add values to the same location twice and add those
#         # dummy_value = ifelse(_LG.src(tr) == _LG.dst(tr),zero(1/dt),oneunit(1/dt)) # dummy_value voegt nu een 1 toe wanneer dest and source gelijk zijn, maar dit moet enkel in geval ni en nj ook nog gelijk zijn. Extra controle in de tijdslus.
#         dst = get_prop(std, tr, :distr) 
#         lb = floor(cquantile(dst, tol) / dt) * dt
#         for (ni,nt) in enumerate(t)
#             if ni > steps + 1
#                 Φ = nt:-dt:t[1]
#                 # Φ represents the sojourn time, ranging from
#                 lt = t[1]:dt:nt
#                 NΦ = length(Φ)
                
#                 append!(I, (Ns * (ni-1) + _LG.dst(tr)).*ones(Int,NΦ))
#                 append!(J, [Ns * (nj-1) + _LG.src(tr) for nj in 1:NΦ])
#                 # append!(V, .- dt .* weights(ni)[1:NΦ] .* (pdf.(dst, nt.-Φ, Φ) .|> unit(dt)^-1))
#                 append!(V, .- dt .* weights(ni)[1:NΦ] .* (pdf.(dst, Φ, lt) .|> unit(dt)^-1))
#             end
#     end end
#     V[isnan.(V)] .= 0.0
#     return _SA.sparse(I, J, V)
# end

# stochastic process
# function solve!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
#     tsim::Number=1.0u"yr", dt::Number=1u"d", tol::Real=1e-8)
#     # get the input
#     dt = dt |> unit(tsim)
#     t = zero(dt):dt:tsim
#     Nt  = length(t)

#     # solve the problem
#     Φ   = zeros(Nt, ns(std))
#     QTI = Quantity{Float64, 𝐓 ^-1, Unitful.FreeUnits{(unit(1/dt)), 𝐓 ^-1, nothing}}[]


#     println("set_U")
#     @time U = set_U(std,t,tol)
#     println("set_A")
#     @time A = ustrip(set_A(std,t))
#     H=U\A*unit(1/t[1])
    
#     # println("set_P")
#     # @time set_P(std, t, H, tol)
#     println("Linear interpolation of H")
#     @time h = [_INT.LinearInterpolation(collect(t), map(x->H[ns(std) * (x-1) + st], 1:Nt)) for st in states(std)]; # splice id H = st:NS:end
#     h_sol = [map(x->H[ns(std) * (x-1) + st], 1:Nt) for st in states(std)];

#     println("Integration of H to Φ")
#     @time for st in states(std)
#         for (ni,nt) in enumerate(t)
#             w   = weights(ni)
#             # TOM: φ <<< t, zero could be higher
#             # l   = zero(dt):dt:nt
#             l = t[1]:dt:nt .|> unit(tsim)
#     #         # NB: ccdf(t-l,φ) where φ = 0.0, GLENN, additional clarification
#     #         # The probability of being in a state st at time nt depends on two things:
#     #         # First, the probability of being in that state initially, times the probability of not having transitioned out of that state until time nt.
#     #         # The probability of not having transitioned until time nt is characterized by the complementary cumulative density function evaluated at the time
#     #         # of entering the state (0) to determine the weight and the sojourn time (nt-0).
#     #         # Second, on the probability of having transitioned to state st at time x and the probability of not having transitioned out of that state during the sojourn time (nt-x).
#     #         # The probability of having transitioned to state st at time x is characterized by the integral of the transition frequency density h of state st.
#     #         # The probability of not having transitioned out of that state during the sojourn time (nt-x) is characterized by the complementary cumulative density function evaluated at 
#     #         # The last transition time x and the sojourn time nt-x.
#             Φ[ni,st] += get_prop(std, st, :init) * ccdf(std, st, nt, zero(dt))
#             Φ[ni,st] += _QGK.quadgk(x -> h[st](x) * ccdf(std, st, nt-x, x), t[1],nt,rtol=1e-7)[1] 
                                
#             # Φ[ni,st] += sum(dt .* w[nj] .* H[ns(std) * (nj-1) + st] .* 
#             #                     ccdf(std, st, nt-nl, nl) 
#             #                     for (nj,nl) in enumerate(l))
#     end end

    
#     # set the output
#     set_prop!(std, :cls, cls)
#     set_prop!(std, :time, t)
#     set_prop!(std, states(std), :prob, [Φ[:,ns] for ns in states(std)])

#     # set the solved status
#     set_info!(std, :solved, true)
#     return h_sol
# end

function solveP!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
    tsim::Number=1.0u"yr", dt::Number=1u"hr", tol::Real=1e-8)
    # get the input
    dt = dt |> unit(tsim)
    t = zero(dt):dt:tsim
    Nt  = length(t)

    # solve the problem
    Φ   = zeros(Nt, ns(std))

    println("set_U")
    @time U = set_U(std,t,tol)
    println("set_A")
    @time A = ustrip(set_A(std,t))
    H=U\A*unit(1/t[1])

    println("set_P")
    @time set_P(std, t, H, tol)
    
    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)

    # set the solved status
    set_info!(std, :solved, true)
end

# function solved!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
#     tsim::Number=1.0u"yr", dt::Number=1u"hr", acc::Int64=10, steps::Number=100.0, tol::Real=1e-8)
    
#     # Get the input for solving numerical inrush phase
#     dt_n = dt/acc |> unit(tsim);
#     t_n = zero(dt_n):dt_n:dt_n*(steps)*acc;

#     U = set_U(std,t_n[1:end],tol)
#     A = ustrip(set_A(std,t_n[1:end]))
#     H1 = U\A*unit(1/t_n[1])
 
#     H_c = compress(H1, steps, std)

#     H1o = [ map(x->H1[ns(std) * (x-1) + st], 1:length(t_n)-ns(std)) for st in states(std)];
#     H1co = [ map(x->H_c[ns(std) * (x-1) + st], 1:Int(length(H_c)/ns(std))) for st in states(std)];

#     # get the input for final problem solution
#     dt = dt |> unit(tsim)
#     t = zero(dt):dt:tsim
#     Nt  = length(t)

#     # solve the problem
#     Φ   = zeros(Nt, ns(std))

#     U = set_U(std,t,steps,tol)
#     A = ustrip(set_A(std,t,H_c,steps))
#     H=U\A*unit(1/t[1])

#     Ho = [ map(x->H[ns(std) * (x-1) + st], 1:Nt) for st in states(std)];

#     # set_P(std, t, H, tol)
#     h = [_INT.LinearInterpolation(collect(t), map(x->H[ns(std) * (x-1) + st], 1:Nt)) for st in states(std)]; # splice id H = st:NS:end

#     for st in states(std)
#         for (ni,nt) in enumerate(t)
#             w   = weights(ni)
#             # TOM: φ <<< t, zero could be higher
#             # l   = zero(dt):dt:nt
#             l = t[1]:dt:nt .|> unit(tsim)
#             # NB: ccdf(t-l,φ) where φ = 0.0, GLENN, additional clarification
#             # The probability of being in a state st at time nt depends on two things:
#             # First, the probability of being in that state initially, times the probability of not having transitioned out of that state until time nt.
#             # The probability of not having transitioned until time nt is characterized by the complementary cumulative density function evaluated at the time
#             # of entering the state (0) to determine the weight and the sojourn time (nt-0).
#             # Second, on the probability of having transitioned to state st at time x and the probability of not having transitioned out of that state during the sojourn time (nt-x).
#             # The probability of having transitioned to state st at time x is characterized by the integral of the transition frequency density h of state st.
#             # The probability of not having transitioned out of that state during the sojourn time (nt-x) is characterized by the complementary cumulative density function evaluated at 
#             # The last transition time x and the sojourn time nt-x.
#             Φ[ni,st] += get_prop(std, st, :init) * ccdf(std, st, nt, zero(dt))
#             Φ[ni,st] += _QGK.quadgk(x -> h[st](x) * ccdf(std, st, nt-x, x), t[1],nt,rtol=1e-7)[1] 
                                
#             # Φ[ni,st] += sum(dt .* w[nj] .* H[ns(std) * (nj-1) + st] .* 
#             #                     ccdf(std, st, nt-nl, nl) 
#             #                     for (nj,nl) in enumerate(l))
#     end end

    
#     # set the output
#     set_prop!(std, :cls, cls)
#     set_prop!(std, :time, t)
#     set_prop!(std, states(std), :prob, [Φ[:,ns] for ns in states(std)])

#     # set the solved status
#     set_info!(std, :solved, true)

#     return H1o, H1co, Ho
# end

function set_P(std::AbstractSTD, t::StepRangeLen, H::Vector, tol::Real)
    Ns  = ns(std)
    Nt = length(t)

    for st in states(std)
        dst_v = [get_prop(std, _LG.Edge(st,nx),:distr) for nx in _LG.outneighbors(std.graph, st)]
        init   = get_prop(std, st, :init)
        h = map(x->H[ns(std) * (x-1) + st], 1:Nt)        
        if init > 0.0
            p = init .*  1-sum(cdf(dstr, 0, zero(t[1])) for dstr in dst_v) .+ integral(dst_v, t, h, tol)
        else
            p = integral(dst_v, t, h, tol)        
        end
        set_prop!(std, st, :prob, p)
    end
end



"""
    MultiStateSystems.weights(x::Int)

Determine integration weights based on extended Simpson's rule. 
w[1] and w[end] = 1/3, even weights = 4/3 and uneven weights = 2/3.
"""
# function weights(x::Int)
#     x==1 && return [0]
#     x==2 && return [1/2, 1/2]
#     x==3 && return [1/3, 4/3, 1/3]
#     x==4 && return [3/8, 9/8, 9/8, 3/8]
#     x==5 && return (2/45) .* [7, 32, 12, 32, 7]
#     x==6 && return (5/288) .* [19, 75, 50, 50, 75, 19]
#     x==7 && return (1/140) .* [41, 216, 27, 272, 27, 216, 41]
#     x==8 && return (7/17280) .* [751, 3577, 1323, 2989, 2989, 1323, 3577, 751]
    
#     weights             = 48 * ones(x)
#     weights[1:4]        = [17, 59, 43, 49]
#     weights[end-3:end]  = [49, 43, 59, 17]
#     return (1/48) .* weights
# end 

# function weights(x::Int)
#     w = zeros(x)
#     for i in 1:length(w)
#         if i % 2 == 0
#            w[i] = 2/3 
#         else
#             w[i] = 4/3
#     end end
#     w[1] = 1/3
#     w[end] = 1/3

#     return w
# end

function weights(N::Int, p::Int)
    if N==1 
        return 0.0
    elseif p == 1 || p == N
        return 0.5
    else
        return 1.0
    end
end

# function weights(x::Int)
#     w = ones(x)
#     return w
# end


# # TOM: Alternative formulation for the state ccdf, using the trapping info prop
# cdf(std::AbstractSTD, ns::Int, φ::Number, t::Number) = 
#     sum(cdf(get_prop(std, _LG.Edge(ns,nx), :distr), φ, t) for nx in _LG.outneighbors(std.graph, ns))
# ccdf(std::AbstractSTD, ns::Int, φ::Number, t::Number)  = 
#     ifelse(get_prop(std, ns, :trapping), 1.0, 1.0 - cdf(std, ns, φ, t))


# TOM: added functionality to enable \ with units.
elunit(A::Matrix{U}) where U = _UF.unit(U)
elunit(b::Vector{U}) where U = _UF.unit(U)
_LA.:\(A::Matrix{<:Number}, b::Vector{<:Number}) = 
    (ustrip.(A) \ ustrip(b)) * (elunit(b) / elunit(A))

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
    # the power sources are repaired at time t. With a lognormal distribution for 
    # the repair time with mean μ and standard deviation σ.
    p_failure = 1 - ntw_av[i]
    p_repair_time_ge_T = ccdf(Exponential(θ),T)
    return 1-p_failure*p_repair_time_ge_T
end

function integral(dst_v::Vector, t::StepRangeLen, h::Vector, tol::Real)
    # controleer voor schaalfactor
    d2h = abs.(diff(diff(ustrip(h))))
    id = 1
    idx = [1]

    while true
        next_id = findfirst(x -> x > tol, cumsum(d2h[id:end]))
        next_id == nothing ? break : ~ ;

        id += next_id
        push!(idx,id)
    end
    push!(idx,length(d2h)+2)
    τ   = t[idx]
    γ   = _INT.interpolate((τ,), h[idx], _INT.Gridded(_INT.Linear()))
    ϕ   = [_QGK.quadgk(x -> γ(x) * (1-sum(cdf(dstr, nτ-x, x) for dstr in dst_v)), t[1], nτ, rtol=tol)[1] for nτ in τ]
    # niet gelijk gespacete interpolatie 
    # Y = interpolate((x,), y, Gridded(Linear()))
    p   = _INT.LinearInterpolation(τ, ϕ)(t)
    println(size(idx))
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