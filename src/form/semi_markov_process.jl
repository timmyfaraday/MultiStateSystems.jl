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

# stochastic process

function solve!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; tsim::Number=1.0u"yr", dt::Number=1.0u"d", tol::Real=1e-8)
    t = zero(dt):dt:tsim
    Nt = length(t)

    # calculate H
    A = zeros(ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        for tr in transitions(std)
            id      = ns(std) * (ni-1) + _LG.dst(tr)
            # met wat moet hier vermeningvuldigd worden?
            # - pure pdf
            # - integral van de pdf over dt
            # - cdf
            A[id]   = get_prop(std, _LG.src(tr), :init) * 
                        pdf(get_prop(std, tr, :distr),0.0,nt)
    end end 
            
    U = zeros(ns(std) * Nt, ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        w = weights(ni)
        # hier kan code optimalisatie gebeuren gegeven dat φ <<< t
        φ = 0.0:dt:nt
        for (nj,nφ) in enumerate(φ)
            Ψ = zeros(ns(std), ns(std))
            for tr in transitions(std)
                Ψ[_LG.dst(tr),_LG.src(tr)] = pdf(get_prop(std, tr, :distr),nφ,nt)
            end
            # ranges for time and sojourn time in U
            rT = (ns(std) * (ni - 1) + 1):(ns(std) * ni)
            rΦ = (ns(std) * (nj - 1) + 1):(ns(std) * nj)
            if ni == nj
                U[rT,rΦ] = Matrix(1.0I, ns(std), ns(std)) .- dt .* w[nj] .* Ψ
            else
                U[rT,rΦ] = -dt .* w[nj] .* Ψ
            end
    end end

    H = U \ A

    Φ = zeros(Nt, ns(std))
    for st in states(std)
        for (ni,nt) in t
            w = weights(ni)
            # hier kan code optimalisatie gebeuren gegeven dat φ <<< t
            φ = 0.0:dt:nt
            Φ[ni,st] += get_prop(std,st,:init) * ccdf(std,st,0.0,nt)
            Φ[ni,st] += sum(dt .* w[nj] .* H[ns(std) * (ni-1) + st] * ccdf(std,st,nφ,nt) for (nj,nφ) in enumerate(φ))
    end end

    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)
    set_prop!(std, states(std), :prob, [Φ[:,ns] for ns in states(std)])
end

#determine integration weights based on extended Simpson's rule. w[1] and w[end] = 1/3, even weights = 4/3 and uneven weights = 2/3.
function weights(x)
    w = zeros(x)
    for i in 1:x
        if i % 2 == 0
           w[i] = 4/3 
        else
            w[i] = 2/3
    end end
    w[1] = 1/3
    w[end] = 1/3

    return w
end

# check whether sum of ccdf could be suffient
ccdf(std::AbstractSTD, ns::Int, φ::Real, t::Real) = 
    1.0 - sum(cdf(get_prop(std,_LG.edge(ns,nx),:distr), φ, t) 
                for nx in _LG.outneighbors(G,ns))





