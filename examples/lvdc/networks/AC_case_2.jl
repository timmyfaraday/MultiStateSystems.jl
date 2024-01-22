################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# load pkgs
using Unitful
using MultiStateSystems

# pkg const
const _MSS = MultiStateSystems

# setting for a specific analysis
cls = SteadyStateProcess()

# include the stochastic models for the elements
stdᵃᶜ   = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/ac.jl"))
stdᵗᶠ   = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/transfo.jl"))
stdⁱⁿᵛ  = include(joinpath(_MSS.BASE_DIR,"examples/lvdc/elements/inv.jl"))


# solve the state-transitions diagrams
solve!(stdᵃᶜ, cls)
solve!(stdᵗᶠ, cls)
solve!(stdⁱⁿᵛ, cls)

# initialize the network corresponding to the low voltage dc system 
ntwˡᵛ = Network()

add_users!(ntwˡᵛ, node = [2,4])
add_sources!(ntwˡᵛ, node = [1], 
                    name = ["ac grid"], 
                    std  = [stdᵃᶜ])
add_components!(ntwˡᵛ, edge = [(1,2),(2,3),(3,4)],
                       name = ["transfo", "ac/dc","dc/ac"],
                       std  = [stdᵗᶠ, stdⁱⁿᵛ, stdⁱⁿᵛ])

solve!(ntwˡᵛ)

ntwˡᵛ.usr[1][:ugf]
ntwˡᵛ.usr[2][:ugf]