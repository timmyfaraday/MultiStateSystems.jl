################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

## Van Acker
# Sets
function set_cycles!(std::AbstractSTD)
    # Initialize a vector of cycles, i.e., vector of edges
    C   = Vector{Vector{_LG.Edge}}()
    # Find the initial state of the problem
    n   = findfirst(x -> x .== 1.0, get_sprop(std,:init)) 
    # Fill the vector of cycles, each beginning in the initial state 
    for nc in _LG.simplecycles(std.graph)
        id = findfirst(x -> x == n, nc)
        ns = push!(circshift(nc, 1-id), n)
        push!(C,[_LG.Edge(ns[ni],ns[ni+1]) for ni in 1:length(ns)-1])
    end
    set_prop!(std, :n, n)
    set_prop!(std, :C, C)
end
function set_failures!(std::AbstractSTD)
    # Set the failure and preventive maintenance decision sets :F and :P
    set_prop!(std, :F, Set([nt for nt in transitions(std) 
                               if get_prop(std, nt, :type) == :f]))
    set_prop!(std, :P, Set([nt for nt in transitions(std) 
                               if get_prop(std, nt, :type) == :p]))
    # Set the failure and preventive maintenance decision subsets :Fᵐⁱⁿ, :Fᵖᵉʳ
    # and :Pᵖᵉʳ
    set_prop!(std, :Fᵐⁱⁿ, Set([nc[1] for nc in get_prop(std,:C)
                                     if get_prop(std, nc[end], :type) == :cmm]))
    set_prop!(std, :Fᵖᵉʳ, Set([nc[1] for nc in get_prop(std,:C)
                                     if get_prop(std, nc[end], :type) == :cpm]))
    set_prop!(std, :Pᵖᵉʳ, Set([nc[1] for nc in get_prop(std,:C)
                                     if get_prop(std, nc[end], :type) == :ppm]))
end
function set_cohorts!(std::AbstractSTD, tsim::Number, tol::Real)
    # Determine the initial sojourn time φinit
    n     = get_prop(std, :n)
    φinit = get_prop(std, n, :φinit)
    # Determine the appropriate cohort step da, cohort set A
    da    = _ST.mean(_UF.maximum(cquantile(get_prop(std, nt, :distr),tol) 
                                            for nt in nc[2:end]) 
                                            for nc in get_prop(std, :C))/2e1
    A     = ceil(Int,-φinit/da):ceil(Int,tsim/da)
    set_prop!(std, :A, A)
    set_prop!(std, :da, da)
    set_prop!(std, :time, zero(da):da:A[end]*da)
    set_prop!(std, n, :φinit, -A[1] * da) 
end
# Parameters
check_minimal(std::AbstractSTD, C::Vector{_LG.Edge}, nf::_LG.Edge) =
    C[1] == nf && get_prop(std, C[end], :type) ∈ [:cmm]
check_perfect(std::AbstractSTD, C::Vector{_LG.Edge}, nf::_LG.Edge) =
    C[1] == nf && get_prop(std, C[end], :type) ∈ [:cpm, :ppm]
function set_sojourn_time!(std::AbstractSTD, tol::Real)
    # Pre-allocate the necessary input
    n     = get_prop(std, :n)
    tsim  = get_prop(std, :time)[end]
    φinit = get_prop(std, n, :φinit)

    # Get an approximate for the overall failure rate
    φmax     = zero(tsim)
    pdf_ovr  = (φ) -> 0.0/_UF.oneunit(φmax)
    ccdf_ovr = (φ) -> 0.0
    for nf in get_prop(std, :F)
        dst  = get_prop(std, nf, :distr)
        φmax = max(φmax, cquantile(dst,tol))
        pdf_ovr  = let  pdf_ovr = pdf_ovr; 
                        (φ) -> pdf_ovr(φ) + pdf.(dst,φ); end
        ccdf_ovr = let  ccdf_ovr = ccdf_ovr;
                        (φ) -> ccdf_ovr(φ) + ccdf.(dst,φ); end
    end
    set_prop!(std, n, :φmax, φmax)
    set_prop!(std, n, :rate, (φ) -> pdf_ovr.(φ)./ccdf_ovr.(φ))
    # Find an appropriate sojourn time series, using an explicit Bogacki-Shampine
    # method of order 3
    φmax  = min(get_prop(std, n, :φmax), φinit + tsim)
    prb   = _ODE.ODEProblem((u,p,t) -> -get_prop(std, n, :rate)(t).*u,
                            1.0, (zero(φmax),φmax))
    φ     = _ODE.solve(prb, _ODE.BS3(), reltol=tol, abstol=tol).t
    # Update the sojourn time series to include the initial sojourn time 
    if !in(φinit,φ)
        id = searchsortedfirst(φ, φmax/5)
        deleteat!(φ,findall(x -> φinit .<= x .<=φinit+φ[id], φ))  
        sort!(push!(φ, (φinit .+ φ[1:id])...)) 
    end
    # Update the sojourn time series to include the preventive maintenance 
    # decision sojourn times.
    for np in get_prop(std, :P)
        φprev = get_prop(std, np, :distr).o
        dφ    = _UF.oneunit(φprev)
        φprev == φmax ? sort!(push!(φ, [φprev-dφ, φprev]...)) : ~ ;
        φprev <  φmax ? sort!(push!(φ, [φprev-dφ, φprev, φprev+dφ]...)) : ~ ;
    end
    # Determine the indices of the preventive maintenance decisions
    ϕprev = Dict(np => findfirst(x -> x == get_prop(std, np, :distr).o,φ) 
                                                for np in get_prop(std, :P))
    set_prop!(std, n, :φ, φ)
    set_prop!(std, n, :φmax, φmax)
    set_prop!(std, n, :ϕprev, ϕprev)
    set_prop!(std, n, :ϕinit, findfirst(x->x.==φinit,φ))
end
function convolute_pmf!(pmf::Vector{Float64}, cycle::Vector{_LG.Edge},
                        std::AbstractSTD, dφ::Number, tol::Real)
    offs, temp, wght = 0.0unit(dφ), ones(Float64,1), one(Float64)
    for nt in filter(x -> !isa(get_prop(std, x, :distr),Dirac), cycle)
        dist = get_prop(std, nt, :distr)
        wght *= weight(dist)
        temp = _DSP.conv(temp,diff(cdf.(dist,sojourn(dist,dφ,tol)))) 
    end
    for nt in filter(x ->  isa(get_prop(std, x, :distr),Dirac), cycle)
        dist = get_prop(std, nt, :distr)
        offs += offset(dist)
        temp *= weight(dist)
        wght *= weight(dist)
    end
    shift!(temp, uconvert(unit(dφ),offs), dφ)
    temp[end] += wght - sum(temp)
    update!(pmf,temp)    
end    
function set_convoluted_pmfs!(std::AbstractSTD, tol::Real)
    da = get_prop(std, :da)
    for nf in get_prop(std, :Fᵐⁱⁿ)
        C = filter(x -> check_minimal(std, x[:], nf), get_prop(std, :C))
        if !isempty(C)
            dφ  = _UF.maximum(_UF.maximum(cquantile(get_prop(std, nt, :distr),tol) 
                                                            for nt in nc[2:end]) 
                                                            for nc in C)/1e6
            pmf = Vector{Float64}()
            for nc in C 
                convolute_pmf!(pmf, nc[2:end], std, dφ, tol)
            end
            φ   = zero(dφ):dφ:length(pmf)*dφ
            cdf = _INT.LinearInterpolation(collect(φ), vcat(0.0,cumsum(pmf)),
                                                extrapolation_bc = _INT.Flat())
            set_prop!(std, nf, :fᵐⁱⁿ, diff(cdf.(zero(da):da:ceil(φ[end]/da)*da)))
    end end
    for nf in get_prop(std, :Fᵖᵉʳ)
        C = filter(x -> check_perfect(std, x[:], nf), get_prop(std, :C))
        if !isempty(C)
            dφ  = _UF.maximum(_UF.maximum(cquantile(get_prop(std, nt, :distr),tol) 
                                                            for nt in nc[2:end]) 
                                                            for nc in C)/1e6
            pmf = Vector{Float64}()
            for nc in C 
                convolute_pmf!(pmf, nc[2:end], std, dφ, tol)
            end
            φ   = zero(dφ):dφ:length(pmf)*dφ
            cdf = _INT.LinearInterpolation(collect(φ), vcat(0.0,cumsum(pmf)), 
                                                extrapolation_bc = _INT.Flat())
            set_prop!(std, nf, :fᵖᵉʳ, diff(cdf.(zero(da):da:ceil(φ[end]/da)*da)))
    end end
    for np in get_prop(std, :Pᵖᵉʳ)
        C = filter(x -> check_perfect(std, x[:], np), get_prop(std, :C))
        if !isempty(C)
            dφ  = _UF.maximum(_UF.maximum(cquantile(get_prop(std, nt, :distr),tol) 
                                                            for nt in nc[2:end]) 
                                                            for nc in C)/1e6
            pmf = Vector{Float64}()
            for nc in C 
                convolute_pmf!(pmf, nc[2:end], std, dφ, tol)
            end
            φ   = zero(dφ):dφ:length(pmf)*dφ
            cdf = _INT.LinearInterpolation(collect(φ), vcat(0.0,cumsum(pmf)), 
                                                extrapolation_bc = _INT.Flat())
            set_prop!(std, np, :fᵖᵉʳ, diff(cdf.(zero(da):da:ceil(φ[end]/da)*da)))
    end end
end
function set_failure_pmfs!(std::AbstractSTD, tol::Real)
    da  = get_prop(std, :da)
    for nf in get_prop(std, :F)
        dst = get_prop(std, nf, :distr)
        φ   = zero(da):da:cquantile(dst,tol)
        set_prop!(std, nf, :f, diff(cdf.(dst,φ)))
    end
end
function set_failure_rates!(na::Int, std::AbstractSTD, tol::Real)  
    # Pre-allocate the necessary input
    n   = get_prop(std, :n)
    da  = get_prop(std, :da)
    𝓕   = get_prop(std, :F)
    tsim = get_prop(std, :time)[end]
    φinit = get_prop(std, n, :φinit)
    φmax = get_prop(std, n, :φmax)

    φx  = ifelse(na<0, ceil(φinit/da)*da, 0.0)
    tx  = ceil(min(φmax,φinit+tsim)/da)*da
    φ   = φx:da:(tx - na*da); Φ = Int((φ[1]+da)/da):Int(φ[end]/da); Nφ = Φ[end]
    t   = (φx + na*da):da:tx; T = Int((t[1]+da)/da):Int(t[end]/da); Nt = length(t)

    # Determine the pmfs of the failures
    f   = Dict(nf => zeros(Float64, Nt-1) for nf in 𝓕)
    for nf in 𝓕 f[nf] += get_prop(std, nf, :f)[Φ] end 
    # Determine the overall survival function R
    R   = ones(Float64, Nt) - [0.0;sum(cumsum(f[nf]) for nf in 𝓕)]
    R   = max.(R, tol.*ones(Float64, Nt))
    # Set the failure rates with respect to da
    for nf in 𝓕 set_prop!(std, nf, :λᵀ, (1/da).*[f[nf][1];f[nf]./R[2:end]]) end
    set_prop!(std, n, :λᵀ, sum(get_prop(std, nf, :λᵀ) for nf in get_prop(std, :F)))
    # Set the failure rates with respect to φ
    # for nf in 𝓕 set_prop!(std, nf, :λᵠ, _INT.LinearInterpolation(collect(t), get_prop(std, nf, :λᵀ))(φ)) end
    # set_prop!(std, n, :λᵠ, sum(get_prop(std, nf, :λᵠ) for nf in get_prop(std, :F)))
    # # Set the interpolation indices
    # i = [find(x ->x >= tol, abs.(diff(get_prop(std, nf, :λᵀ)))) for nf in 𝓕]
    # for np in get_prop(std, :P) push!(i,floor(Int,get_prop!(std, nr, distr).o/da)) end
    # set_prop!(std, n, :int, i)
end