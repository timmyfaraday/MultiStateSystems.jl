################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
"""
Example considering a multi-state Markov model for a coal fired power generating
unit taken from:

> A Multi-State Markov Model for a Short-Term Reliability Analysis of a Power
  Generating Unit by A. Lisnianski, D. Elmakias, and H. Ben Haim (2012)

A unit's available generating capacity can take any real value in [0,gⁿᵒᵐ].
This is approximated by a discrete state continuous time stochastic process
Gₐ(t) as follows:
- Two states to describe the limits of the generating capacity:
    - State 1, where g₁ = 0
    - State N, where gₙ = gⁿᵒᵐ
- N-2 additional states seperated with an 'averaged' output gᵢ within an
  interval Δg = gⁿᵒᵐ/(N-2).

If the random times between transitions are ignored and only the transition time
instances are of interest, the stochastic process Gₐ(t) reduces to a Markov
model Gₐₘ(t), with approximated transition intensities âᵢⱼ from state i to
state j≠i.

The case study in the paper considers a four-state model to represent a coal
fired power generating unit, with the input given in Table 1[^1].

Table 1: States' output and transition intensities.
| g₁  = 0 MW        | g₂  = 247 MW      | g₃  = 482 MW      | g₄  = 576 MW      |
|-------------------|-------------------|-------------------|-------------------|
| â₁₂ ≈ 0.08000 1/h | â₂₁ ≈ 0.02912 1/h | â₃₁ ≈ 0.00000 1/h | â₄₁ ≈ 0.00015 1/h |
| â₁₃ ≈ 0.01333 1/h | â₂₃ ≈ 0.32353 1/h | â₃₂ ≈ 0.02885 1/h | â₄₂ ≈ 0.00010 1/h |
| â₁₄ ≈ 0.00000 1/h | â₂₄ ≈ 0.02912 1/h | â₃₄ ≈ 0.35577 1/h | â₄₃ ≈ 0.00069 1/h |

[^1]: The data in the table is an approximation. Exact results are obtained by
using the fraction of kᵢⱼ over Tᵢ as given in the paper.
"""

# load pkgs
using Unitful
using MultiStateSystems

# initialize the state-transition diagram corresponding to the coal fired power
# generating unit.
stdᶜᶠ = STD()

# add the states to the std
add_states!(stdᶜᶠ, power = [0u"MW", 247u"MW", 482u"MW", 576u"MW"],
                   init  = [0.0, 0.0, 0.0, 1.0])

# add the transitions to the std
add_transitions!(stdᶜᶠ, states = [(1,2),(1,3)],
                        rate = [0.08000u"1/hr",0.01333u"1/hr"])
add_transitions!(stdᶜᶠ, states = [(2,1),(2,3),(2,4)],
                        rate = [0.02941u"1/hr",0.32353u"1/hr",0.02941u"1/hr"])
add_transitions!(stdᶜᶠ, states = [(3,2),(3,4)],
                        rate = [0.02885u"1/hr",0.35577u"1/hr"])
add_transitions!(stdᶜᶠ, states = [(4,1),(4,2),(4,3)],
                        rate = [0.00015u"1/hr",0.00010u"1/hr",0.00069u"1/hr"])

# return the std
return stdᶜᶠ