################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
# Contributor: Gayan Abeynayake ([@gayan86](https://github.com/gayan86))       #
################################################################################

function cluster_wind_power(input::Vector{Float64}; number_of_clusters::Int=10)
    if any(x->x.<=0.0,input) input[input.<=0.0] .= 0.0 end

    temp = copy(input)
    number_of_samples = length(input)
    clusters = Jenks.JenksClassification(number_of_clusters-2,input,errornorm=2)

    output = round.([minimum(temp),clusters.centres...,maximum(temp)],digits=3)
    bounds = clusters.bounds
    assign = zeros(Int,number_of_samples)
    assign[temp.==bounds[1]] .= 1
    for nc in 2:number_of_clusters-1
        assign[bounds[nc-1].<temp.<=bounds[nc]] .= nc
    end
    assign[temp.==bounds[end]] .= number_of_clusters

    rate = zeros(Float64,number_of_clusters,number_of_clusters)
    for ni in 1:number_of_samples-1 rate[assign[ni],assign[ni+1]] += 1 end
    rate ./= (sum(rate,dims=2)/number_of_samples)
    rate[_LA.diagind(rate)] .= 0

    return output, rate
end
