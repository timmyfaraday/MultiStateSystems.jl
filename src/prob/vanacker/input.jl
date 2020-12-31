################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

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
    ns    = findfirst(x -> x .== 1.0, get_sprop(std,:init)) 
    φinit = get_prop(std, ns, :φinit)
    # Determine the appropriate cohort step da, cohort set A
    da    = _ST.mean(maximum(cquantile(get_prop(std, nt, :distr),tol) 
                                            for nt in nc[2:end]) 
                                            for nc in get_prop(std, :C))/2e2
    A     = ceil(Int,-φinit/da):ceil(Int,tsim/da)
    set_prop!(std, :A, A)
    set_prop!(std, :da, da)
end

# # Parameters
# check_minimal(std::AbstractSTD, C::Vector{_LG.Edge}, nf::_LG.Edge) =
#     C[1] == nf && get_prop(std, C[end], :type) ∈ [:cmm]
# check_perfect(std::AbstractSTD, C::Vector{_LG.Edge}, nf::_LG.Edge) =
#     C[1] == nf && get_prop(std, C[end], :type) ∈ [:cpm, :ppm]
# function set_sojourn_time!(std::AbstractSTD, tsim::Number, tol::Real)
#     # Pre-allocate the necessary input
#     n     = get_prop(std, :n)
#     φinit = get_prop(std, n, :φinit)

#     # Get an approximate for the overall failure rate
#     φmax     = zero(tsim)
#     pdf_ovr  = (φ) -> 0.0/_UF.oneunit(φmax)
#     ccdf_ovr = (φ) -> 0.0
#     for nf in get_prop(std, :F)
#         dst  = get_prop(std, nf, :distr)
#         φmax = max(φmax, cquantile(dst,tol))
#         pdf_ovr  = let  pdf_ovr = pdf_ovr; 
#                         (φ) -> pdf_ovr(φ) + pdf.(dst,φ); end
#         ccdf_ovr = let  ccdf_ovr = ccdf_ovr;
#                         (φ) -> ccdf_ovr(φ) + ccdf.(dst,φ); end
#     end
#     set_prop!(std, n, :φmax, φmax)
#     set_prop!(std, n, :rate, (φ) -> pdf_ovr.(φ)./ccdf_ovr.(φ))
#     # Find an appropriate sojourn time series, using an explicit Bogacki-Shampine
#     # method of order 3
#     φmax  = min(get_prop(std, n, :φmax), φinit + tsim)
#     prb   = _ODE.ODEProblem((u,p,t) -> -get_prop(std, s, :rate)(t).*u,
#                             1.0, (zero(φmax),φmax))
#     φ     = _ODE.solve(prb, _ODE.BS3(), reltol=tol, abstol=tol).t
#     # Update the sojourn time series to include the initial sojourn time 
#     if !in(φinit,φ)
#         id = searchsortedfirst(φ, φmax/5)
#         deleteat!(φ,findall(x -> φinit .<= x .<=φinit+φ[id], φ))  
#         sort!(push!(φ, (φinit .+ φ[1:id])...)) 
#     end
#     # Update the sojourn time series to include the preventive maintenance 
#     # decision sojourn times.
#     for np in get_prop(std, :P)
#         φprev = get_prop(std, np, :distr).o
#         dφ    = _UF.oneunit(φprev)
#         φprev == φmax ? sort!(push!(φ, [φprev-dφ, φprev]...)) : ~ ;
#         φprev <  φmax ? sort!(push!(φ, [φprev-dφ, φprev, φprev+dφ]...)) : ~ ;
#     end
#     # Determine the indices of the preventive maintenance decisions
#     ϕprev = Dict{_LG.Edge,Int}(np => findfirst(x -> x == get_prop(std, np, :distr).offset,φ) 
#                                                     for np in get_prop(std, :P))
#     set_prop!(std, ns, :φ, φ)
#     set_prop!(std, ns, :φmax, φmax)
#     set_prop!(std, ns, :ϕprev, ϕprev)
#     set_prop!(std, ns, :ϕinit, findfirst(x->x.==φinit,φ))
# end
# function set_rates!(na::Int, std::AbstractSTD, tol::Real)
#     ## Initialize the necessary parameters
#     # Settings parameters
#     tᵐᵃˣ = s.tᵐᵃˣ
#     dt   = s.dt
#     εⁿ   = s.εⁿ
#     # Feeder parameters
#     L    = ϕ.attr["𝓛"]
#     Na   = ϕ.attr["Na"]
#     Lⁿ   = ϕ.attr["𝓛ⁿ"]
#     Lᵃ   = ϕ.attr["𝓛ᵃ"]
#     Lᶠ   = ϕ.attr["𝓛ᶠ"]
#     φⁱⁿⁱ = ϕ.attr["φⁱⁿⁱ"]
#     φᵐᵃˣ = ϕ.attr["φᵐᵃˣ"]
#     φᵖʳᵉ = ϕ.attr["φᵖʳᵉ"]
#     # STD parameters
#     Tr   = ϕ.std.trans
#     δ    = ϕ.std.attr["δ"]
#     F    = ϕ.std.attr["𝓕"]
#     Fʷ   = ϕ.std.attr["𝓕ʷ"]
#     Rᵖ   = ϕ.std.attr["𝓡ᵖ"]
    
#     ## Initialize the necessary parameters
#     # Maximum index of the failure set, and the maximum index of the age set
#     Nf   = maximum(F);
#     # Determine the initial normal sojourn time xⁱⁿⁱ and intial time zⁱⁿⁱ
#     xⁱⁿⁱ = ifelse(na<0,ceil(φⁱⁿⁱ/dt)*dt,0.0)
#     zⁱⁿⁱ = xⁱⁿⁱ + na*dt
#     # Determine the maximal normal sojourn time xᵐᵃˣ and time zᵐᵃˣ
#     xᵐᵃˣ = ceil(min(φᵐᵃˣ,φⁱⁿⁱ+tᵐᵃˣ)/dt)*dt - na*dt
#     zᵐᵃˣ = xᵐᵃˣ + na*dt
#     # Initialize the sojourn time φ, age α and time t and their respective ranges
#     φ    = xⁱⁿⁱ:dt:xᵐᵃˣ; Φ = Int((φ[1]+dt)/dt):Int((φ[end])/dt); Nφ = Φ[end]
#     t    = zⁱⁿⁱ:dt:zᵐᵃˣ; T = Int((t[1]+dt)/dt):Int((t[end])/dt); Nt = length(t)
    
#     # Pre-allocate output
#     f   = Dict(nf => zeros(Float64, length(t)) for nf in get_prop(std, :F))
    
#     for nf in get_prop(std, :F) 
#     if !δ[0].α
#         # Set the pmf of all non-load dependent failures
#         for nf ∈ F if nf ∉ Fʷ f[nf] += Tr[nf].attr["f"][Φ] end end
#         # Set the pmf of all load dependent failures
#         for nf ∈ F if nf ∈ Fʷ for nl ∈ Lᶠ[nf], nα ∈ keys(L[nl].attr["fʷ"])
#             nl ∈ Lⁿ ? f[nf] += L[nl].attr["p"][nα][T].*L[nl].attr["fʷ"][nα][Φ] : ~ ;
#             nl ∈ Lᵃ ? f[nf] += L[nl].attr["p̂"][nα][T].*L[nl].attr["fʷ"][nα][Φ] : ~ ;
#         end end end
#     else
#         # Initialize the survival function R
#         R  = 1.0
#         # Initialize a from- fr, to- to and value array vl as well as the initial
#         # state probability vector p
#         fr = []; to = [1]; vl = []; p  = zeros(Float64,Na+Nf); p[Φ[1]] = 1.0
#         # Enumerate over all sojourn time indices and set the pmf of all load-dependent
#         # failures
#         for nφ in Φ[1:end-1]
#             # Determine the Markov matrix at sojourn time index nφ
#             fr, to, vl = get_Transitions(na,nφ,Na,R,unique(to),p,s,ϕ)
#             P          = sparse(to,fr,vl,Na+Nf,Na+Nf)
#             # Determine the state probability vector p at sojourn time index nφ and
#             # update the failure pmf f
#             p          = P*p
#             for nf ∈ F if nf ∈ F f[nf][nφ] = p[Na+nf] end end
#             # Update the survival function
#             for nf ∈ F if nf ∈ F R -= p[Na+nf] end end
#         end
#     end
#     # Determine the discrete failure rates Λ and associated functions λ
#     i = Array{Int,1}()
#     λ = Dict{Int,Function}()
#     Λ = Dict{Int,Array{Float64,1}}(); for nf ∈ F Λ[nf] = zeros(Float64,Nt+na) end
#     # Determine the overall survival function R
#     R = ones(Float64,Nt) - [0.0;sum(cumsum(f[nf]) for nf in F)]
#     R = max.(R,εⁿ.*ones(Float64,Nt))
#     # Determine the individual discrete failure rates (dt) and failure functions
#     for nf ∈ F Λ[nf][1:Nt] += (1/dt).*[f[nf][1];f[nf]./R[2:end]] end
#     for nf ∈ F λ[nf] = φn -> Spline1D(φ,Λ[nf][1:Nt];k=3,bc="nearest")(φn) end
#     # Determine the overall discrete failure rate (dt) and failure function
#     Λ[0] = sum(Λ[nf] for nf ∈ F)
#     λ[0] = φn -> Spline1D(φ,Λ[0];k=3,bc="nearest")(φn)
#     # Determine the indices at which additional interpolation is necessary
#     for nf ∈ F vcat(i,find(x->x>=εⁿ,abs.(diff(Λ[nf])))) end
#     for nr ∈ Rᵖ push!(i,floor(Int,φᵖʳᵉ[nr]/dt)) end
    
#     # Set the failure rate dictionaries of the feeder ϕ
#     ϕ.std.attr["λ"]  = λ
#     ϕ.std.attr["λⁱ"] = i
#     ϕ.std.attr["λᵀ"] = Λ
    
#     end
# end
# function convolute_pmf!(pmf::Vector{Float64}, cycle::Vector{_LG.Edge},
#                         std::AbstractSTD, dφ::Number, tol::Real)
#     offs, temp, wght = 0.0unit(dφ), ones(Float64,1), one(Float64)
#     for nt in filter(x -> !isa(get_prop(std, x, :distr),Dirac), cycle)
#         dist = get_prop(std, nt, :distr)
#         wght *= weight(dist)
#         temp = _DSP.conv(temp,diff(cdf.(dist,sojourn(dist,dφ,tol)))) 
#     end
#     for nt in filter(x ->  isa(get_prop(std, x, :distr),Dirac), cycle)
#         dist = get_prop(std, nt, :distr)
#         offs += offset(dist)
#         temp *= weight(dist)
#         wght *= weight(dist)
#     end
#     shift!(temp, uconvert(unit(dφ),offs), dφ)
#     temp[end] += wght - sum(temp)
#     update!(pmf,temp)    
# end    
# function set_convoluted_pmfs!(std::AbstractSTD, tol::Real)
#     da = get_prop(std, :da)
#     for nf in get_prop(std, :Fᵐⁱⁿ)
#         C = filter(x -> check_minimal(std, x[:], nf), get_prop(std, :C))
#         if !isempty(C)
#             dφ  = maximum(maximum(cquantile(get_prop(std, nt, :distr),tol) 
#                                                         for nt in nc[2:end]) 
#                                                         for nc in C)/1e6
#             pmf = Vector{Float64}()
#             for nc in C 
#                 convolute_pmf!(pmf, nc[2:end], std, dφ, tol)
#             end
#             φ   = (0.0)unit(dφ):dφ:length(pmf)*dφ
#             cdf = _INT.LinearInterpolation(collect(φ), vcat(0.0,cumsum(pmf)))
#             φᵃ  = (0.0)unit(da):da:ceil(φ[end]/da)*da
#             pmf = cdf(φᵃ[2:end])-cdf(φ)
#             set_prop!(std, nf, :φᵐⁱⁿ, φ[end])
#             set_prop!(std, nf, :fᵐⁱⁿ, pmf)
#     end end
#     for nf in get_prop(std, :Fᵖᵉʳ)
#         C = filter(x -> check_perfect(std, x[:], nf), get_prop(std, :C))
#         if !isempty(C)
#             dφ  = maximum(maximum(cquantile(get_prop(std, nt, :distr),tol) 
#                                                         for nt in nc[2:end]) 
#                                                         for nc in C)/1e6
#             pmf = Vector{Float64}()
#             for nc in C 
#                 convolute_pmf!(pmf, nc[2:end], std, dφ, tol)
#             end
#             cdf  = vcat(0.0,cumsum(pmf))
#             φ    = (0.0)unit(dφ):dφ:length(pmf)*dφ
#             itp  = _INT.LinearInterpolation(collect(φ), cdf)
            
#             set_prop!(std, nf, :φᵖᵉʳ, φ[end])
#             set_prop!(std, nf, :fᵖᵉʳ, _INT.LinearInterpolation(collect(φ), cdf))
#     end end
# end

# # function set_rates!(na::Int, std::AbstractSTD, tol::Real; trans=transitions(std))
# #     ## Initialize the necessary parameters
# #     # Settings parameters
# #     tᵐᵃˣ = s.tᵐᵃˣ
# #     dt   = s.dt
# #     εⁿ   = s.εⁿ
# #     # Feeder parameters
# #     L    = ϕ.attr["𝓛"]
# #     Na   = ϕ.attr["Na"]
# #     Lⁿ   = ϕ.attr["𝓛ⁿ"]
# #     Lᵃ   = ϕ.attr["𝓛ᵃ"]
# #     Lᶠ   = ϕ.attr["𝓛ᶠ"]
# #     φⁱⁿⁱ = ϕ.attr["φⁱⁿⁱ"]
# #     φᵐᵃˣ = ϕ.attr["φᵐᵃˣ"]
# #     φᵖʳᵉ = ϕ.attr["φᵖʳᵉ"]
# #     # STD parameters
# #     Tr   = ϕ.std.trans
# #     δ    = ϕ.std.attr["δ"]
# #     F    = ϕ.std.attr["𝓕"]
# #     Fʷ   = ϕ.std.attr["𝓕ʷ"]
# #     Rᵖ   = ϕ.std.attr["𝓡ᵖ"]
    
# #     ## Initialize the necessary parameters
# #     # Maximum index of the failure set, and the maximum index of the age set
# #     Nf   = maximum(F);
# #     # Determine the initial normal sojourn time xⁱⁿⁱ and intial time zⁱⁿⁱ
# #     xⁱⁿⁱ = ifelse(na<0,ceil(φⁱⁿⁱ/dt)*dt,0.0)
# #     zⁱⁿⁱ = xⁱⁿⁱ + na*dt
# #     # Determine the maximal normal sojourn time xᵐᵃˣ and time zᵐᵃˣ
# #     xᵐᵃˣ = ceil(min(φᵐᵃˣ,φⁱⁿⁱ+tᵐᵃˣ)/dt)*dt - na*dt
# #     zᵐᵃˣ = xᵐᵃˣ + na*dt
# #     # Initialize the sojourn time φ, age α and time t and their respective ranges
# #     φ    = xⁱⁿⁱ:dt:xᵐᵃˣ; Φ = Int((φ[1]+dt)/dt):Int((φ[end])/dt); Nφ = Φ[end]
# #     t    = zⁱⁿⁱ:dt:zᵐᵃˣ; T = Int((t[1]+dt)/dt):Int((t[end])/dt); Nt = length(t)
    
# #     ## MAIN
# #     # Determine pmfs (dt) associated with all failures nf ∈ 𝓕 along the cohort na,
# #     # depending on whether the overall failure rate is degradation dependent
# #     f = Dict{Int,Array{Float64,1}}(); for nf ∈ F f[nf] = zeros(Float64,Nφ) end
# #     if !δ[0].α
# #         # Set the pmf of all non-load dependent failures
# #         for nf ∈ F if nf ∉ Fʷ f[nf] += Tr[nf].attr["f"][Φ] end end
# #         # Set the pmf of all load dependent failures
# #         for nf ∈ F if nf ∈ Fʷ for nl ∈ Lᶠ[nf], nα ∈ keys(L[nl].attr["fʷ"])
# #             nl ∈ Lⁿ ? f[nf] += L[nl].attr["p"][nα][T].*L[nl].attr["fʷ"][nα][Φ] : ~ ;
# #             nl ∈ Lᵃ ? f[nf] += L[nl].attr["p̂"][nα][T].*L[nl].attr["fʷ"][nα][Φ] : ~ ;
# #         end end end
# #     else
# #         # Initialize the survival function R
# #         R  = 1.0
# #         # Initialize a from- fr, to- to and value array vl as well as the initial
# #         # state probability vector p
# #         fr = []; to = [1]; vl = []; p  = zeros(Float64,Na+Nf); p[Φ[1]] = 1.0
# #         # Enumerate over all sojourn time indices and set the pmf of all load-dependent
# #         # failures
# #         for nφ in Φ[1:end-1]
# #             # Determine the Markov matrix at sojourn time index nφ
# #             fr, to, vl = get_Transitions(na,nφ,Na,R,unique(to),p,s,ϕ)
# #             P          = sparse(to,fr,vl,Na+Nf,Na+Nf)
# #             # Determine the state probability vector p at sojourn time index nφ and
# #             # update the failure pmf f
# #             p          = P*p
# #             for nf ∈ F if nf ∈ F f[nf][nφ] = p[Na+nf] end end
# #             # Update the survival function
# #             for nf ∈ F if nf ∈ F R -= p[Na+nf] end end
# #         end
# #     end
# #     # Determine the discrete failure rates Λ and associated functions λ
# #     i = Array{Int,1}()
# #     λ = Dict{Int,Function}()
# #     Λ = Dict{Int,Array{Float64,1}}(); for nf ∈ F Λ[nf] = zeros(Float64,Nt+na) end
# #     # Determine the overall survival function R
# #     R = ones(Float64,Nt) - [0.0;sum(cumsum(f[nf]) for nf in F)]
# #     R = max.(R,εⁿ.*ones(Float64,Nt))
# #     # Determine the individual discrete failure rates (dt) and failure functions
# #     for nf ∈ F Λ[nf][1:Nt] += (1/dt).*[f[nf][1];f[nf]./R[2:end]] end
# #     for nf ∈ F λ[nf] = φn -> Spline1D(φ,Λ[nf][1:Nt];k=3,bc="nearest")(φn) end
# #     # Determine the overall discrete failure rate (dt) and failure function
# #     Λ[0] = sum(Λ[nf] for nf ∈ F)
# #     λ[0] = φn -> Spline1D(φ,Λ[0];k=3,bc="nearest")(φn)
# #     # Determine the indices at which additional interpolation is necessary
# #     for nf ∈ F vcat(i,find(x->x>=εⁿ,abs.(diff(Λ[nf])))) end
# #     for nr ∈ Rᵖ push!(i,floor(Int,φᵖʳᵉ[nr]/dt)) end
    
# #     # Set the failure rate dictionaries of the feeder ϕ
# #     ϕ.std.attr["λ"]  = λ
# #     ϕ.std.attr["λⁱ"] = i
# #     ϕ.std.attr["λᵀ"] = Λ
    
# #     end
    