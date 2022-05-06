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

# parameters
function set_parameters!(std::AbstractSTD, cls::SemiMarkovProcess)
    
end

function set_A(std::AbstractSTD, t::StepRange, dt::Number, Nt::Int64)
    "A indicates that state j can be reached if the process was initially
    in state i and remains there until time t." 

    A = zeros(ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        for tr in transitions(std)
            id      = ns(std) * (ni-1) + _LG.dst(tr)
            A[id]   = A[id] + get_prop(std, _LG.src(tr), :init)unit(dt) * 
                        pdf(get_prop(std, tr, :distr), nt - 0.0u"hr",0.0u"hr")
    end end 
    replace_nan!(A)
    return A
end

function set_U(std::AbstractSTD, t::StepRange, dt::Number, Nt::Int64)
    "U indicates that state j can be reached if the process enters state i 
    at time φ and remains there until time t"
    U = zeros(ns(std) * Nt, ns(std) * Nt)

    for (ni,nt) in enumerate(t)
        w = weights(ni)
        # hier kan code optimalisatie gebeuren gegeven dat φ <<< t
        φ = zero(dt):dt:nt
        for (nj,nφ) in enumerate(φ)
            Ψ = zeros(ns(std), ns(std))unit(1/dt)
            for tr in transitions(std)
                Ψ[_LG.dst(tr),_LG.src(tr)] = pdf(get_prop(std, tr, :distr),nt-nφ,nφ)
                # distr betreden op nϕ en blijven tot nt. Sojourn time = nt-nϕ
            end
            # ranges for time and state entrance time in U
            rT = (ns(std) * (ni - 1) + 1):(ns(std) * ni)
            rΦ = (ns(std) * (nj - 1) + 1):(ns(std) * nj)
            if ni == nj
                U[rT,rΦ] = Matrix(1.0_LA.I, ns(std), ns(std)) .- dt .* w[nj] .* Ψ
            else
                U[rT,rΦ] = -dt .* w[nj] .* Ψ
            end
    end end

    replace_nan!(U)
    return U
end

# stochastic process
function solve!(std::AbstractSTD, cls::AbstractSemiMarkovProcess; 
                tsim::Number=1.0u"yr", dt::Number=1.0u"d", tol::Real=1e-8)
    # set the input
    set_parameters!(std,cls)    

    # get the input
    t = zero(dt):dt:tsim
    Nt = length(t)

    # solve the problem
    # calculate H

    H = set_U(std, t, dt, Nt)\set_A(std, t, dt, Nt)

    Φ = zeros(Nt, ns(std))
    for st in states(std)
        for (ni,nt) in enumerate(t)
            w = weights(ni)
            # hier kan code optimalisatie gebeuren gegeven dat φ <<< t
            φ = zero(dt):dt:nt
            Φ[ni,st] += get_prop(std,st,:init) * ccdf(std, st, nt - 0.0u"hr",0.0u"hr")
            Φ[ni,st] += sum(dt .* w[nj]unit(1/dt) .* H[ns(std) * (nj-1) + st] * ccdf(std, st, nt - nφ, nφ) for (nj,nφ) in enumerate(φ))
    end end

    # set the output
    set_prop!(std, :cls, cls)
    set_prop!(std, :time, t)
    set_prop!(std, states(std), :prob, [Φ[:,ns] for ns in states(std)])

    # set the solved status
    set_info!(std, :solved, true)
end

#determine integration weights based on extended Simpson's rule. w[1] and w[end] = 1/3, even weights = 4/3 and uneven weights = 2/3.
function weights(x)
    # w = zeros(x)
    # for i in 1:x
    #     if i % 2 == 0
    #        w[i] = 4/3 
    #     else
    #         w[i] = 2/3
    # end end
    # w[1] = 1/3
    # w[end] = 1/3

    # return w
    if x==1
        weights = [0]
    elseif x==2
        weights = [1/2;1/2]
    elseif x==3
        weights = [1/3;4/3;1/3]
    elseif x==4
        weights = [3/8;9/8;9/8;3/8]
    elseif x==5
        weights = (2/45).*[7;32;12;32;7]
    elseif x==6
        weights = (5/288).*[19;75;50;50;75;19]
    elseif x==7 
        weights = (1/140).*[41;216;27;272;27;216;41]
    elseif x==8
        weights = (7/17280).*[751;3577;1323;2989;2989;1323;3577;751]
    else
        weights = 48*ones(x)
        weights[1:4] = [17;59;43;49]
        weights[end-3:end] = [49;43;59;17]
        weights = (1/48).*weights
    end
end 

# check whether sum of ccdf could be suffient
ccdf(std::AbstractSTD, ns::Int, φ::Quantity, t::Quantity)  = 
    1.0 - mysum(std, ns, φ, t)


function mysum(std::AbstractSTD, ns::Int, φ::Quantity, t::Quantity) 
# Summation function to overcome summation error over empty collection.
    if _LG.outneighbors(std.graph, ns) != Int64[]
        F = sum(cdf(get_prop(std,_LG.Edge(ns,nx),:distr), φ, t)
            for nx in _LG.outneighbors(std.graph, ns))
    else
        F = 0
end end

function replace_nan!(x)
    for i = eachindex(x)
        if isnan(x[i])
            x[i] = 0
    end end
end


