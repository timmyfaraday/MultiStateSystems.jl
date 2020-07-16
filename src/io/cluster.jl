################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
# Contributor: Gayan Abeynayake (@gayan86)                                     #
################################################################################

"""
# Wind Power
"""
function cluster_wind_power(input::Array; number_of_clusters::Int=10)
    if !isa(input,Matrix{Float64}) input = convert(Matrix{Float64},input) end
    if size(input)[2] < number_of_clusters input = transpose(input) end

    number_of_samples = length(input)
    clusters = _CL.kmeans(input,number_of_clusters; maxiter=200, tol=1e-10)

    output = round.(vec(clusters.centers), digits = 3)
    assign = clusters.assignments

    temp = zeros(Float64,number_of_clusters,number_of_clusters)
    for ni in 1:number_of_samples-1 temp[assign[ni],assign[ni+1]] += 1 end
    temp ./= (sum(temp,dims=2)/number_of_samples)
    temp[_LA.diagind(temp)] .= 0

    idx = sortperm(output)
    rate = zeros(Float64,number_of_clusters,number_of_clusters)
    for nc in CartesianIndices(temp) rate[idx[nc[1]],idx[nc[2]]] = temp[nc] end

    return output[idx], rate
end
