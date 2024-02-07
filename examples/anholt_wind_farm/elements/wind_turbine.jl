################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# cluster function
function cluster_wind_power(input::Matrix{Float64}; number_of_clusters::Int=10)
    if any(x->x.<=0.0,input) input[input.<=0.0] .= 0.0 end

    number_of_samples = length(input)
    temp = vec(reshape(deepcopy(input), number_of_samples, 1))
    clusters = Clustering.kmeans(input, number_of_clusters-2)

    idx = sortperm(vec(clusters.centers))

    output = round.([minimum(temp),clusters.centers[idx]...,maximum(temp)],digits=3)
    assign = [idx[nx] for nx in clusters.assignments] .+ 1
    assign[temp.==0.0] .= 1
    assign[temp.==3.6] .= number_of_clusters

    rate = zeros(Float64,number_of_clusters,number_of_clusters)
    for ni in 1:number_of_samples-1 rate[assign[ni],assign[ni+1]] += 1 end
    rate ./= (sum(rate, dims=2) / number_of_samples)
    rate[LinearAlgebra.diagind(rate)] .= 0
    rate[isnan.(rate)] .= 0.0

    return output, rate
end

## STD
# Read in the wind power output data, and determine the clustered outputs and
# corresponding rates.
path = joinpath(BASE_DIR,"examples/anholt_wind_farm/data/Anholt.xlsx")
wind_power = XLSX.readdata(path,"Sheet1!D2:D52215")
wind_power = convert(Matrix{Float64}, reshape(wind_power, 1, length(wind_power)))
output, rate = cluster_wind_power(wind_power, number_of_clusters = number_of_clusters)
init = zeros(length(output)); init[1] = 1.0

# Initialize the state-transition diagrams corresponding to the output (wto) and
# reliability (wtr) of a wind turbine.
stdʷᵗᵒ = STD()
stdʷᵗʳ = STD()

# Add the states to the std's
add_states!(stdʷᵗᵒ, power = (output)u"MW",
                    init  = init)
mss ? add_states!(stdʷᵗʳ, power = [(Inf)u"MW",0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW",
                                   0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW",0.0u"MW"],
                          init = [1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]) :
      add_states!(stdʷᵗʳ, power = [(Inf)u"MW",0.0u"MW"],
                          init = [1.0,0.0]) ;

# Add the transitions to the std's
add_transitions!(stdʷᵗᵒ, rate = (rate)u"1/yr")
if mss
    add_transitions!(stdʷᵗʳ, states = [(1,2),(2,1)],
                             rate = [0.0590u"1/yr",0.0132u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,3),(3,1)],
                             rate = [0.0070u"1/yr",0.1695u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,4),(4,1)],
                             rate = [0.0770u"1/yr",0.0158u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,5),(5,1)],
                             rate = [0.0420u"1/yr",0.0361u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,6),(6,1)],
                             rate = [0.0240u"1/yr",0.3704u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,7),(7,1)],
                             rate = [0.3380u"1/yr",0.0442u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,8),(8,1)],
                             rate = [0.4320u"1/yr",0.0752u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,9),(9,1)],
                             rate = [0.4370u"1/yr",0.0625u"1/hr"])
    add_transitions!(stdʷᵗʳ, states = [(1,10),(10,1)],
                             rate = [0.5380u"1/yr",0.0515u"1/hr"])
else
    add_transitions!(stdʷᵗʳ, states = [(1,2),(2,1)],
                             rate = [1.954u"1/yr",0.04715u"1/hr"])
end
# Solve the problems
solve!(stdʷᵗᵒ, SteadyStateProcess())
solve!(stdʷᵗʳ, SteadyStateProcess())

## Network
# Initialize the wind turbine network
ntwʷᵗ = Network()

# Add the user, sources and components to the wind turbine network ntwʷᵗᵒ
add_user!(ntwʷᵗ, node = 1)
add_source!(ntwʷᵗ, node = 2, std = stdʷᵗᵒ, dep = true)
wtr ? add_component!(ntwʷᵗ, edge = (2,1), std = stdʷᵗʳ) :
      add_component!(ntwʷᵗ, edge = (2,1)) ;
