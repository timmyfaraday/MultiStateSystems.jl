################################################################################
#  Copyright 2022, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

@testset "Distributions" begin
    Ω(t) = max(0.0,min(1.0,0.01u"1/hr" * t))
    
    @testset "Exponential" begin
        # input 
        θ, ω, φ = 10.0u"hr", 0.8, 3.0u"hr"

        # constructors
        @test 𝑬() == Exponential() == Exponential(1.0,1.0)
        @test 𝑬(θ) == Exponential(θ) == Exponential(θ,1.0)
        @test 𝑬(θ,ω) == Exponential(θ,ω)
        @test 𝑬(θ,Ω) == Exponential(θ,Ω)

        # general functions 
        dst = Exponential(θ,ω)

        @test _MSS.scale(dst) == θ
        @test _MSS.weight(dst) == ω
        @test _MSS.params(dst) == (θ, ω)  
        
        @test _MSS.rate(dst) == 1 / θ

        @test _MSS.minimum(dst) == zero(θ)
        @test _MSS.maximum(dst) == (Inf)unit(θ)

        # quantile functions
        # https://keisan.casio.com/menu/system/000000000540 
        @test isapprox(_MSS.quantile(dst,0.4),5.1082562376u"hr",atol=1e-10u"hr")
        @test _MSS.quantile(dst,0.8) == _MSS.cquantile(dst,0.2)

        # density functions
        # https://keisan.casio.com/menu/system/000000000540 
        @test isapprox(_MSS.pdf(dst,φ),ω * 0.0740818220u"1/hr",atol=1e-10u"1/hr")
        @test isapprox(_MSS.cdf(dst,φ),ω * 0.2591817793,atol=1e-10)
        @test _MSS.cdf(dst,φ) == ω - _MSS.ccdf(dst,φ)
    end

    @testset "Weibull" begin
        # dummy input 
        θ, α, ω, φ = 10.0u"hr", 1.3, 0.8, 3.0u"hr"

        # constructors
        @test 𝑾() == Weibull() == Weibull(1.0,1.0,1.0)
        @test 𝑾(θ) == Weibull(θ) == Weibull(θ,1.0,1.0)
        @test 𝑾(θ,α) == Weibull(θ,α) == Weibull(θ,α,1.0)
        @test 𝑾(θ,α,ω) == Weibull(θ,α,ω)
        @test 𝑾(θ,α,Ω) == Weibull(θ,α,Ω)

        # general functions 
        dst = Weibull(θ,α,ω)

        @test _MSS.scale(dst) == θ
        @test _MSS.shape(dst) == α
        @test _MSS.weight(dst) == ω
        @test _MSS.params(dst) == (θ, α, ω)  

        @test _MSS.minimum(dst) == zero(θ)
        @test _MSS.maximum(dst) == (Inf)unit(θ)

        # quantile functions
        # https://keisan.casio.com/menu/system/000000000540 
        @test isapprox(_MSS.quantile(dst,0.4),5.9647790998u"hr",atol=1e-10u"hr")
        @test _MSS.quantile(dst,0.8) == _MSS.cquantile(dst,0.2)

        # density functions
        # https://keisan.casio.com/menu/system/000000000540 
        @test isapprox(_MSS.pdf(dst,φ),ω * 0.0735002655u"1/hr",atol=1e-10u"1/hr")
        @test isapprox(_MSS.cdf(dst,φ),ω * 0.1886482463,atol=1e-10)
        @test _MSS.cdf(dst,φ) == ω - _MSS.ccdf(dst,φ)
    end

    @testset "Normal" begin
        # dummy input 
        μ, σ, ω, φ = 10.0u"hr", 2.0u"hr", 0.8, 12.0u"hr"

        # constructors
        @test 𝑵() == Normal() == Normal(0.0,1.0,1.0)
        @test 𝑵(μ) == Normal(μ) == Normal(μ,oneunit(μ),1.0)
        @test 𝑵(μ,σ) == Normal(μ,σ) == Normal(μ,σ,1.0)
        @test 𝑵(μ,σ,ω) == Normal(μ,σ,ω)
        @test 𝑵(μ,σ,Ω) == Normal(μ,σ,Ω)

        # general functions 
        dst = Normal(μ,σ,ω)

        @test _MSS.mean(dst) == μ
        @test _MSS.std(dst) == σ
        @test _MSS.var(dst) == σ^2
        @test _MSS.weight(dst) == ω
        @test _MSS.params(dst) == (μ, σ, ω)  

        @test _MSS.minimum(dst) == -Inf * unit(μ)
        @test _MSS.maximum(dst) == Inf * unit(μ)
        @test _MSS.median(dst) == μ

        # quantile functions
        # For Normal(10, 2), quantile at p=0.6: 10 + 2*invnormcdf(0.6) ≈ 10 + 2*0.2533 ≈ 10.5066
        @test isapprox(_MSS.quantile(dst,0.6), 10.5066u"hr", atol=1e-3u"hr")
        @test _MSS.quantile(dst,0.8) == _MSS.cquantile(dst,0.2)
        @test _MSS.quantile(dst,0.5) == _MSS.median(dst)

        # density functions
        # For Normal(10, 2) at x=12: z=(12-10)/2=1, pdf = (1/(2*√(2π))) * exp(-1²/2) ≈ 0.1209854
        @test isapprox(_MSS.pdf(dst,φ), ω * 0.1209854u"1/hr", atol=1e-6u"1/hr")
        # For Normal(10, 2) at x=12: z=1, cdf = (1 + erf(1/√2))/2 ≈ 0.8413447 
        @test isapprox(_MSS.cdf(dst,φ), ω * 0.8413447, atol=1e-6)
        @test _MSS.cdf(dst,φ) == ω - _MSS.ccdf(dst,φ)
    end

    @testset "LogNormal" begin
        # dummy input 
        μ, σ, ω, φ = 1.0u"hr", 0.25u"hr", 0.8, 3.0u"hr"

        # constructors
        @test 𝑳() == LogNormal() == LogNormal(1.0,1.0,1.0)
        @test 𝑳(μ) == LogNormal(μ) == LogNormal(μ,oneunit(μ),1.0)
        @test 𝑳(μ,σ) == LogNormal(μ,σ) == LogNormal(μ,σ,1.0)
        @test 𝑳(μ,σ,ω) == LogNormal(μ,σ,ω)
        @test 𝑳(μ,σ,Ω) == LogNormal(μ,σ,Ω)

        # general functions 
        dst = LogNormal(μ,σ,ω)

        @test _MSS.weight(dst) == ω
        @test _MSS.params(dst) == (μ, σ, ω)  

        @test _MSS.minimum(dst) == zero(μ)
        @test _MSS.maximum(dst) == (Inf)unit(μ)

        # quantile functions
        # https://keisan.casio.com/menu/system/000000000540 
        @test isapprox(_MSS.quantile(dst,0.4),2.5514535699u"hr",atol=1e-10u"hr")
        @test _MSS.quantile(dst,0.8) == _MSS.cquantile(dst,0.2)

        # density functions
        # https://keisan.casio.com/menu/system/000000000540 
        @test isapprox(_MSS.pdf(dst,φ),ω * 0.4921107291u"1/hr",atol=1e-10u"1/hr")
        @test isapprox(_MSS.cdf(dst,φ),ω * 0.6533752704,atol=1e-10)
        @test _MSS.cdf(dst,φ) == ω - _MSS.ccdf(dst,φ)
    end

end