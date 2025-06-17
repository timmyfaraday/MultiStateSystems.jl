################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

# using pkgs
using JSON
using MultiStateSystems

const _MSS = MultiStateSystems

# path to and read-in json
json_file = joinpath(_MSS.BASE_DIR,"examples/lvdc/dcide/data/test.json")
json_data = JSON.parsefile(json_file)

# solve the problem
solve!(json_data)