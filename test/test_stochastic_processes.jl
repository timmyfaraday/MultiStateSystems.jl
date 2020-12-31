# ################################################################################
# #  Copyright 2020, Tom Van Acker                                               #
# ################################################################################
# # MultiStateSystems.jl                                                         #
# # A Julia package to solve multi-state system models.                          #
# # See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
# ################################################################################
# @testset "Stochastic Processes" begin

#     @testset "Markovian System"
#         # Initialize the state-transition diagram
#         stdᵐʳᵏ = STD()
#         # Add the states to the state-transition diagram
#         add_states!(stdᵐʳᵏ, init = [1.0, 0.0, 0.0, 0.0])
#         # Add the transitions to the state-transition diagram
#         add_transitions!(stdᵐʳᵏ, states = [(1,2),(2,3),(2,4),(4,3),(3,1)],
#                                  type   = [:f, :r, :r, :r, :pcm],
#                                  distr  = [𝑬(1.0e+4u"hr",1.0),
#                                            𝑬(1.0e+1u"hr",0.8),
#                                            𝑬(1.0e+1u"hr",0.2),
#                                            𝑬(1.0e+2u"hr",1.0),
#                                            𝑬(1.0e+3u"hr",1.0)])

#         # Solve the state-transition diagram using the different stochastic
#         #processes
#         solᵐʳᵏ = solve(stdᵐʳᵏ,1.0u"yr", alg = :markov_process)
#         solˢᵐᵖ = solve(stdᵐʳᵏ,1.0u"yr", alg = :semimarkov_process)
#         solᵛᵃⁿ = solve(stdᵐʳᵏ,1.0u"yr", alg = :vanacker_process)

#         # Tests
#         for ns in 1:ns(stdᵐʳᵏ)
#             @test  isapprox(solᵐʳᵏ[ns],solˢᵐᵖ[ns]; atol = 1e-6)
#             @test  isapprox(solˢᵐᵖ[ns],solᵛᵃⁿ[ns]; atol = 1e-6)
#             @test  isapprox(solᵛᵃⁿ[ns],solᵐʳᵏ[ns]; atol = 1e-6)
#         end
#     end

#     @testset "Competing Risks"
#         # Initialize the state-transition diagram
#         stdᶜᵖᵗ = STD()
#         # Add the states to the state-transition diagram
#         add_states!(stdᶜᵖᵗ, init = [1.0, 0.0, 0.0, 0.0,
#                                     0.0, 0.0, 0.0,
#                                     0.0, 0.0, 0.0])
#         # Add the transitions to the state-transition diagram
#         add_transitions!(stdᶜᵖᵗ, states = [(1,2),(2,3),(2,4),(4,3),(3,1),
#                                            (1,5),(5,6),(5,7),(7,6),(6,1),
#                                            (1,8),(8,9),(8,10),(10,9),(9,1)],
#                                  type   = [:f, :r, :r, :r, :pcm,
#                                            :f, :r, :r, :r, :pcm,
#                                            :f, :r, :r, :r, :pcm],
#                                  distr  = [𝑬(1.0e+4u"hr",0.3),
#                                            𝑬(1.0e+1u"hr",0.8),
#                                            𝑬(1.0e+1u"hr",0.2),
#                                            𝑬(1.0e+2u"hr",1.0),
#                                            𝑬(1.0e+3u"hr",1.0),
#                                            𝑬(1.0e+4u"hr",0.4),
#                                            𝑬(1.0e+1u"hr",0.8),
#                                            𝑬(1.0e+1u"hr",0.2),
#                                            𝑬(1.0e+2u"hr",1.0),
#                                            𝑬(1.0e+3u"hr",1.0),
#                                            𝑬(1.0e+4u"hr",0.3),
#                                            𝑬(1.0e+1u"hr",0.8),
#                                            𝑬(1.0e+1u"hr",0.2),
#                                            𝑬(1.0e+2u"hr",1.0),
#                                            𝑬(1.0e+3u"hr",1.0)])

#         # Solve the state-transition diagram using the different stochastic
#         # processes
#         solᵐʳᵏ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :markov_process)
#         solˢᵐᵖ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :semimarkov_process)
#         solᵛᵃⁿ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :vanacker_process)

#         # Tests
#         for ns in 1:ns(stdᵐʳᵏ)
#             @test !isapprox(solᵐʳᵏ[ns],solˢᵐᵖ[ns]; atol = 1e-6)
#             @test  isapprox(solˢᵐᵖ[ns],solᵛᵃⁿ[ns]; atol = 1e-6)
#             @test !isapprox(solᵛᵃⁿ[ns],solᵐʳᵏ[ns]; atol = 1e-6)
#         end
#     end

#     @testset "Semi-Markovian System"
#         # Initialize the state-transition diagram
#         stdˢᵐˢ = STD()
#         # Add the states to the state-transition diagram
#         add_states!(stdˢᵐˢ, init = [1.0, 0.0, 0.0, 0.0,
#                                     0.0, 0.0, 0.0,
#                                     0.0, 0.0, 0.0])
#         # Add the transitions to the state-transition diagram
#         add_transitions!(stdˢᵐˢ, states = [(1,2),(2,3),(2,4),(4,3),(3,1),
#                                            (1,5),(5,6),(5,7),(7,6),(6,1),
#                                            (1,8),(8,9),(8,10),(10,9),(9,1)],
#                                  type   = [:f, :r, :r, :r, :pcm,
#                                            :f, :r, :r, :r, :pcm,
#                                            :f, :r, :r, :r, :pcm],
#                                  distr  = [𝑬(0.4e+4u"hr",0.3),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.8),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.2),
#                                            𝑳𝑵(1.0e+2u"hr",0.08u"hr",1.0),
#                                            𝑾(1.0e+3u"hr",10.0,1.0),
#                                            𝑬(1.0e+4u"hr",0.08,0.4),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.8),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.2),
#                                            𝑳𝑵(1.0e+2u"hr",0.08u"hr",1.0),
#                                            𝑾(1.0e+3u"hr",10.0,1.0),
#                                            𝑬(2.5e+4u"hr",0.3),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.8),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.2),
#                                            𝑳𝑵(1.0e+2u"hr",0.08u"hr",1.0),
#                                            𝑾(1.0e+3u"hr",10.0,1.0)])

#         # Solve the state-transition diagram using the different stochastic
#         # processes
#         solᵐʳᵏ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :markov_process)
#         solˢᵐᵖ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :semimarkov_process)
#         solᵛᵃⁿ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :vanacker_process)

#         # Tests
#         for ns in 1:ns(stdᵐʳᵏ)
#             @test !isapprox(solᵐʳᵏ[ns],solˢᵐᵖ[ns]; atol = 1e-6)
#             @test  isapprox(solˢᵐᵖ[ns],solᵛᵃⁿ[ns]; atol = 1e-6)
#             @test !isapprox(solᵛᵃⁿ[ns],solᵐʳᵏ[ns]; atol = 1e-6)
#         end
#     end

#     @testset "Non-Renewal System"
#         # Initialize the state-transition diagram
#         stdˢᵐˢ = STD()
#         # Add the states to the state-transition diagram
#         add_states!(stdˢᵐˢ, init = [1.0, 0.0, 0.0, 0.0,
#                                     0.0, 0.0, 0.0,
#                                     0.0, 0.0, 0.0])
#         # Add the transitions to the state-transition diagram
#         add_transitions!(stdˢᵐˢ, states = [(1,2),(2,3),(2,4),(4,3),(3,1),
#                                            (1,5),(5,6),(5,7),(7,6),(6,1),
#                                            (1,8),(8,9),(8,10),(10,9),(9,1)],
#                                  type   = [:f, :r, :r, :r, :mcm,
#                                            :f, :r, :r, :r, :mcm,
#                                            :f, :r, :r, :r, :pcm],
#                                  distr  = [𝑬(0.4e+4u"hr",0.3),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.8),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.2),
#                                            𝑳𝑵(1.0e+2u"hr",0.08u"hr",1.0),
#                                            𝑾(1.0e+3u"hr",10.0,1.0),
#                                            𝑬(1.0e+4u"hr",0.08,0.4),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.8),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.2),
#                                            𝑳𝑵(1.0e+2u"hr",0.08u"hr",1.0),
#                                            𝑾(1.0e+3u"hr",10.0,1.0),
#                                            𝑬(2.5e+4u"hr",0.3),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.8),
#                                            𝑳𝑵(1.0e+1u"hr",0.08u"hr",0.2),
#                                            𝑳𝑵(1.0e+2u"hr",0.08u"hr",1.0),
#                                            𝑾(1.0e+3u"hr",10.0,1.0)])

#         # Solve the state-transition diagram using the different stochastic
#         # processes
#         solᵐʳᵏ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :markov_process)
#         solˢᵐᵖ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :semimarkov_process)
#         solᵛᵃⁿ = solve(stdᶜᵖᵗ,1.0u"yr", alg = :vanacker_process)

#         # Tests
#         for ns in 1:ns(stdᵐʳᵏ)
#             @test !isapprox(solᵐʳᵏ[ns],solˢᵐᵖ[ns]; atol = 1e-6)
#             @test !isapprox(solˢᵐᵖ[ns],solᵛᵃⁿ[ns]; atol = 1e-6)
#             @test !isapprox(solᵛᵃⁿ[ns],solᵐʳᵏ[ns]; atol = 1e-6)
#         end
#     end
    
# end
