################################################################################
#  Copyright 2025, Tom Van Acker, Glenn Emmers                                 #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################

using MultiStateSystems
using Unitful
using JLD2
using SMTPClient

const _MSS = MultiStateSystems

include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/DC_Faults.jl"))
include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/std_setup.jl"))
include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/lvdc_netw_funcs.jl"))
include(joinpath(_MSS.BASE_DIR,"examples/lvdc/ICDCM_2025/functions/helper.jl"))

# Grid characteristics:
V_DC    = 750; # V
I_max   = 200; # A
V_min   = V_DC*0.85; # V
I²t     = 5000; #A²t
L_p     = 0.0; # H
C_b     = 5.0e-2;# F ->  single bus capacitance
λᶜ = 0.0000743u"1/yr/m"; # Cable failure rate
P_zone = [0.99, 0.999, 0.9999, 1.0]; # Probability of clearing a fault in the first zone
# P_zone = 1.0

# Integration parameters
t_max   = 0.1u"s"
n       = 100000

# Define the characteristics of the protection devices
λ = Dict("SSCB" => [0.05, "CB"], 
         "HCB"  => [0.05, "CB"],
         "MCCB" => [0.05, "CB"],
         "Fuse" => [0.05, "Fuse"])
         
μ = Dict("SSCB" => 11e-6, 
         "HCB"  => 1.1e-3,
         "MCCB" => 5.1e-3,
         "Fuse" => 0)

# Define the outgoing feeder lengths
L_c = Dict("C1" => 1u"m":1u"m":100u"m", 
           "C2" => 1u"m":1u"m":200u"m",
           "C3" => 1u"m":1u"m":150u"m",
           "C4" => 1u"m":1u"m":50u"m",
           "C5" => 1u"m":1u"m":100u"m",)

L_s = Dict("S1" => 1u"m":1u"m":20u"m", 
           "S2" => 1u"m":1u"m":20u"m")

L_b = Dict("C" => 1u"m":1u"m":200u"m")

# Define the protection device used for each feeder
L_load = Dict("SSCB" => L_c,
             "MCCB" => L_c,
             "HCB"  => L_c,
             "Fuse" => L_c,
             "Fuse_MCCB" => L_c,)

L_source = Dict("SSCB" => L_s)

L_bridge = Dict("SSCB" => L_b)

P_cables_loads = calculate_P(L_load, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
P_cables_sources = calculate_P(L_source, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
P_cables_bridge = calculate_P(L_bridge, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)

P_converter_loads = calculate_Pc(L_load, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)
P_converter_sources = calculate_Pc(L_source, L_p, C_b, V_DC, I_max, V_min, n, t_max, μ, λ)

type = "Markov" # "SemiMarkov" or "Markov"
tsim = 25.0u"yr"    # Simulation time
dt = 0.5u"hr"       # Time step

if type == "SemiMarkov"
    cls = SemiMarkovProcess()
elseif type == "Markov"
    cls = MarkovProcess()
else
    error("Invalid type: $type. Use 'SemiMarkov' or 'Markov'.")
end

std_loads = Dict()
std_sources = Dict()
std_bridge = Dict()

for P_z in P_zone
    std_cables_loads = fill_std(L_load, P_cables_loads, λᶜ; type = type, Pc = P_converter_loads, P_zone = P_z)
    std_cables_sources = fill_std(L_source, P_cables_sources, λᶜ; type = type, Pc = P_converter_sources, P_zone = P_z)
    std_cables_bridge = fill_std(L_bridge, P_cables_bridge, λᶜ; type = type, P_zone = P_z)

    solve_std_dict!(std_cables_loads, cls, tsim, dt)
    solve_std_dict!(std_cables_sources, cls, tsim, dt)
    solve_std_dict!(std_cables_bridge, cls, tsim, dt)

    if type == "SemiMarkov"
        std_loads[P_z] = create_reduced_std!(std_cables_loads)
        std_sources[P_z] = create_reduced_std!(std_cables_sources)
        std_bridge[P_z] = create_reduced_std!(std_cables_bridge)
    elseif type == "Markov"
        std_loads[P_z] = std_cables_loads
        std_sources[P_z] = std_cables_sources
        std_bridge[P_z] = std_cables_bridge
    else
        error("Invalid type: $type. Use 'SemiMarkov' or 'Markov'.")
    end

end

send_email_notification()

jldsave(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/std_loads_1bus_Markov.jld"); dict=std_loads, overwrite = true)
jldsave(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/std_sources_1bus_Markov.jld"); dict=std_sources, overwrite = true)
jldsave(joinpath(_MSS.BASE_DIR, "examples/lvdc/ICDCM_2025/results/std_bridge_1bus_Markov.jld"); dict=std_bridge, overwrite = true)

