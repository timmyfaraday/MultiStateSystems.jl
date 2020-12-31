################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# Matrix A
function set_A(na::Int, std::AbstractSTD)
    # Pre-allocate the necessary input
    n   = get_prop(std, :n)
    φ   = get_prop(std, n, :φ)
    λ   = get_prop(std, n, :rate_dφ)
    na > 0 ? ni = 1 : ni = get_prop(std, n, :ϕinit) ;

    # Constructing the sparse matrix A for the trapezoidal integration method
    A   = spdiagm( (vcat(ones(ni),(diff(φ[ni:end]).*λ[ni+1:end]+2)./2),         # First diagonal [1,(h(φ)*f(φ)+2)/2] : n elements
                    vcat(zeros(ni-1),(diff(φ[ni:end]).*λ[ni:end-1]-2)./2)),     # Second diagonal [(h(φ)*f(φ-1)-2)/2] : n-1 elements
                    [0,-1])        
end
function update_A!(na::Int, A::SparseMatrixCSC{Float64,Int64}, std::AbstractSTD)
    # Pre-allocate the necessary input
    n   = get_prop(std, :n)
    φ   = get_prop(std, n, :φ)
    λ   = get_prop(std, n, :rate_dφ)
    na > 0 ? ni = 1 : ni = get_prop(std, n, :ϕinit) ;

    # Constructing the sparse matrix A for the trapezoidal integration method
    A[:] = spdiagm((vcat(ones(ni),(diff(φ[ni:end]).*λ[ni+1:end]+2)./2),         # First diagonal [1,(h(φ)*f(φ)+2)/2] : n elements
                    vcat(zeros(ni-1),(diff(φ[ni:end]).*λ[ni:end-1]-2)./2)),     # Second diagonal [(h(φ)*f(φ-1)-2)/2] : n-1 elements
                    [0,-1])   
end

# Vector b 
function set_b(na::Int, i::Float64, nh::Array{Float64,1}, std::AbstractSTD)
    # Pre-allocate the necessary input
    n   = get_prop(std, :n)
    φ   = get_prop(std, n, :φ)
    na > 0 ? ni = 1 : ni = get_prop(std, n, :ϕinit)
    # Pre-allocate the output
    b   = zeros(Float64, length(φ))

    # Combine the intial value i and non-homogeneous contribution nh
    b[ni]       = i
    b[ni+1:end] .= diff(φ[ni:end]) .* (nh[ni:end-1] + nh[ni+1:end]) ./2         # Midpoint method: h(t)*(f(t-1)+f(t))/2
end
function update_b!(na::Int, i::Float64, nh::Array{Float64,1}, 
                   b::Array{Float64,1}, std::AbstractSTD)
    # Pre-allocate the necessary input
    n   = get_prop(std, :n)
    φ   = get_prop(std, n, :φ)
    na > 0 ? ni = 1 : ni = get_prop(std, n, :ϕinit)

    # Combining the initial value and the non-homogeneous contribution
    b[ni]       = i                                                         
    b[ni+1:end] .= diff(φ[ni:end]) .* (nh[ni:end-1] + nh[ni+1:end]) ./2         # Midpoint method: h(t)*(f(t-1)+f(t))/2
end

# Vector p
function adjust_p!(A::SparseMatrixCSC{Float64,Int64}, pa::Array{Float64,1},
                   qa::Array{Float64,1}, pm::Dict{Int,Float64}, std::AbstractSTD)
    # Pre-allocate the necessary input
    𝓟       = get_prop(std, :P)
    n       = get_prop(std, :n)
    ϕᵖʳᵉᵛ   = get_prop(std, n, :ϕprev)
    # Pre-allocate temporary variables
    b       = zeros(Float64,length(pa))

    qa[:] = pa[:]
    for nr in 𝓟
        ni = ϕᵖʳᵉᵛ[nr] 
        pm[nr] = pa[ni]*get_prop(std, nr, :distr).ω 
        # Update the qa-array
        b[ni] = pm[nr]
        qa    .-= A\b
        # Reset the b-array
        b[ni] = 0.0
    end
end

function interpolate_p!(na::Int, pa::Array{Float64,1}, Pa::Array{Float64,1}
                        std::AbstractSTD)
    # Pre-allocate the necessary input
    da  = get_prop(std, :da)
    n   = get_prop(std, :n)
    φ   = get_prop(std, n, :φ) 
    na > 0 ? ni = 1 : ni = get_prop(std, n, :ϕinit)

    # 
    φᵃ  = φ[ni]:da:φ[end]
    np  = length(Pa)-length(φᵃ)
    Pa[1:np]     = 0.0
    Pa[np+1:end] = evaluate(Spline1D(φ[ni:end],pa[ni:end];k=3),φᵃ)
end

function update_pf!(na::Int, pa::Array{Float64,1}, Pa::Array{Float64,1},
                    pf::Dict{Int,Array{Float64,1}}, std::AbstractSTD)
    # Pre-allocate the necessary input
    𝓕   = get_prop(std, :F)
    da  = get_prop(std, :da)
    n   = get_prop(std, :n)
    φinit = get_prop(std, n, :φinit)

    λᶲ   = ϕ.std.attr["λᶲ"] %%% 
    λᵀ   = ϕ.std.attr["λᵀ"] %%%
    λⁱ   = ϕ.std.attr["λⁱ"] %%%

    na >= 0 ? nt = 1 : nt = ceil(Int,φinit/da) ;
    na >= 0 ? nϕ = 1 : nϕ = get_prop(std, n, :ϕinit) ;
    Pf   = zeros(Float64,length(Pa))
    
    for nf in 𝓕
        Pf = λᵀ[nf].*Pa
        pf[nf][nt:end] = da./2.0.*(Pf[nt:end-1]+Pf[nt+1:end])
        if !isempty(λⁱ)
            Pᶲ = Spline1D(Φ[nϕ:end],λᶲ[nf][nϕ:end].*pa[nϕ:end];k=3)
            for ni in λⁱ
                fr, to = (ni-1)*da, ni*da
                fr >= φinit ? pf[nf][ni] = integrate(Pᶲ,fr,to) : ~ ;
            end
        end
    end
end

function update_P!(na::Int, Pa::Array{Float64,1}, Pn::Array{Float64,1},
                   pf::Dict{Int,Array{Float64,1}}, Pf::Dict{Int,Array{Float64,1}},
                   pp::Dict{Int,Float64}, Pp::Dict{Int,Array{Float64,1}},
                   tsim::Number, std::AbstractSTD)
    # Pre-allocate the necessary input
    da      = get_prop(std, da)
    n       = get_prop(std, :n)
    φmax    = get_prop(std, n, :φmax)
    φinit   = get_prop(std, n, :φinit)

    # Update the initial state probability Pn
    φmax    = ceil(min(φmax, φinit + tsim) / da) * da
    RΦ      = na * da:da:φmax + na * da
    RP      = findfirst(x -> x.>=(0.0)unit(da), RΦ):findlast(x -> x.<=tsim, RΦ)
    RPn     = max(1, na+1):max(1, na+1) + length(RP) - 1
    Pn[RPn] += Pa[RP]
    # Update the failure probability dictionary Pf
    for nf in get_prop(std, :F)
        Pf[nf][RPf] += pf[nf][RP]
    end
    # Update the preventive maintenance dictionary Pp
    for np in get_prop(std, :P)
        ni = ceil(Int, na + (get_prop(std, np, :distr).ω / da))
        if ni < length(Pp[nr]) Pp[nr][ni] += pp[nr] end
    end
end

OA = OffsetArrays.OffsetArray{Float64,1,Array{Float64,1}}

function update_i!(na::Int, pf::Dict{Int,Array{Float64,1}}, pp::Dict{Int,Float64},
                   i::OA, std::AbstractSTD)
    # Pre-allocate the necessary input
    da      = get_prop(std, :da)

    # Update the initial values vector
    for nf in get_prop(std, :Fᵖᵉʳ)
        f     = get_prop(std, nf, :fᵖᵉʳ)
        cnv   = vcat(0.0,convn(pf[nf],f))
        Ri    = na+1:min(linearindices(i)[end],na+length(cnv))
        Rcnv  = 1:min(length(cnv),linearindices(i)[end]-na)
        i[Ri] += cnv[Rcnv]
    end
    for nr in get_prop(std, :Pᵖᵉʳ)
        f     = get_prop(std, nf, :fᵖᵉʳ)
        ni    = ceil(Int,na+φᵖʳᵉ[nr]/da) %%%%%
        Ri    = ni:ni+length(f)-1
        Rp    = 1:min(length(f),linearindices(i)[end]-ni+1)
        i[Ri] += pp[nr].*f
    end
end

# Matrix nh
function reset_nh!(nh::Array{Float64,2})
    # Rotate the nh-matrix with one position, and reset its last vector
    nh        = circshift(nh,(0,-1))
    nh[:,end] = zeros(size(nh,1))
end
function update_nh!(pa::Array{Float64,1}, nh::Array{Float64,2}, std)
    # Update the non-homogeneous matrix
    for nf in get_prop(std, :Fᵐⁱⁿ)
        f   = get_prop(std, nf, :fᵐⁱⁿ)
        nh[:,1:size(NH,2)] += (λᶲ[nf].*p)*f[nf]' %%%%
    end
end