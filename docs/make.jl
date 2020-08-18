################################################################################
#  Copyright 2020, Tom Van Acker                                               #
################################################################################
# MultiStateSystems.jl                                                         #
# A Julia package to solve multi-state system models.                          #
# See http://github.com/timmyfaraday/MultiStateSystems.jl                      #
################################################################################
using Documenter
using MultiStateSystems

makedocs(
    modules     = [MultiStateSystems],
    format      = Documenter.HTML(analytics = "UA-367975-10", mathengine = Documenter.MathJax()),
    sitename    = "MultiStateSystems.jl",
    authors     = "Tom Van Acker",
    pages       = [ "Home"              => "index.md",
                    "Getting Started"   => "quickguide.md",
                    "DSL Manual"        =>
                        [ "State-Transition Diagram"=> "std.md",
                          "Distributions"           => "distribution.md",
                          "Generating Functions"    => "ugf.md",
                          "Network"                 => "network.md"],
                    "Models Manual"     =>
                        [ "Stochastic Processes"    => "processes.md"
                          "Generating Operators"    => "ugo.md"]
                  ]
)

deploydocs(
    repo = "github.com/timmyfaraday/MultiStateSystems.jl.git"
)
