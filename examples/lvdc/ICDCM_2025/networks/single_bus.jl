using MultiStateSystems
using Unitful
using JLD2

const _MSS = MultiStateSystems

jld = jldopen(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/1bus_Markov.jld"), "r")

std_loads = jld["loads"]
std_sources = jld["sources"]
std_bridge = jld["bridge"]

