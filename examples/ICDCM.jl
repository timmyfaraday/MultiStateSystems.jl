using MultiStateSystems
using Unitful
using Plots

const _MSS = MultiStateSystems
P_timely = 0.95;
P_zone   = 0.99;


tsim = 25.0u"yr"  # Simulation time
dt = 0.1u"hr"     # Time step

t = 0.0u"yr":dt:tsim   # Time vector in years

cls = MarkovProcess()

std_feeder = STD()
add_states!(std_feeder, name  = ["A", "U1", "U2","U3", "V2","V3"],
                        power = [(Inf)u"MW", 0.0u"MW", (Inf)u"MW", (Inf)u"MW", 0.0u"MW", 0.0u"MW"],
                        init  = [1.0, 0.0, 0.0, 0.0, 0.0, 0.0])
add_transitions!(std_feeder, states = [(1,2),(1,3),(1,4),(2,1),(3,5),(4,6),(5,1),(6,1)],
                        distr = [Exponential(20.0u"yr", P_timely*P_zone),
                                 Exponential(20.0u"yr", (1-P_timely)*P_zone),
                                 Exponential(20.0u"yr", 1-P_zone),
                                 Exponential(10.0u"d"),
                                 Exponential(10.0u"hr"),
                                 Exponential(10.0u"hr"),
                                 Exponential(10.0u"d"),
                                 Exponential(10.0u"d")
                                ])

solve!(std_feeder, cls; tsim, dt)


ntw_test = Network()

stdˢ = solvedSTD(prob = [ones(length(t))], time = collect(t), power = [(Inf)u"MW"]);
add_source!(ntw_test, node = 1, std = stdˢ);
add_components!(ntw_test, edge = [(1,2), (2,3)],
                        std = [std_feeder, std_feeder]);

Uᵦ = sum(_MSS.get_sprop(nf, :prob)[ns] for nf in [std_feeder, std_feeder] for ns in 3:4)
stdᵇ = solvedSTD(prob = [1 .- Uᵦ, Uᵦ], power = [(Inf)u"MW", 0.0u"MW"],
                 time = collect(t))
add_component!(ntw_test, node = 2, std = stdᵇ)
add_user!(ntw_test, node = 3)

solve!(ntw_test)

Zone_1 = [std_feeder, std_feeder]
Zone_2 = [std_feeder, std_feeder]

ntw_test_2 = Network()
add_source!(ntw_test_2, node = 1, std = stdˢ)
add_source!(ntw_test_2, node = 2, std = stdˢ)
add_components!(ntw_test_2, edge = [(1,3), (2,4), (3,4),(4,3), (3,5), (4,6)],
                        std = [std_feeder, std_feeder, stdˢ, stdˢ, std_feeder, std_feeder])

Uᵦ₁ = sum(_MSS.get_sprop(nf, :prob)[ns] for nf in Zone_1 for ns in 3:4)
stdᵇ¹ = solvedSTD(prob = [1 .- Uᵦ₁, Uᵦ₁], power = [(Inf)u"MW", 0.0u"MW"], time = collect(t))
Uᵦ₂ = sum(_MSS.get_sprop(nf, :prob)[ns] for nf in Zone_2 for ns in 3:4)
stdᵇ² = solvedSTD(prob = [1 .- Uᵦ₂, Uᵦ₂], power = [(Inf)u"MW", 0.0u"MW"], time = collect(t))

add_component!(ntw_test_2, node = 3, std = stdᵇ¹)
add_component!(ntw_test_2, node = 4, std = stdᵇ²)
add_users!(ntw_test_2, node = [3,4,5,6])

solve!(ntw_test_2)

ntw_test_2.usr[1][:ugf].prb[2][end]
ntw_test_2.usr[2][:ugf].prb[2][end]
ntw_test_2.usr[3][:ugf].prb[2][end]
ntw_test_2.usr[4][:ugf].prb[2][end]
