using Documenter
using MultiStateSystems

makedocs(
    modules     = [MultiStateSystems],
    format      = Documenter.HTML(mathengine = Documenter.MathJax()),
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
                        [ "Dependence"              => "dependence.md"
                          "Stochastic Processes"    => "processes.md"
                          "Generating Operators"    => "ugo.md"]
                  ]
)

deploydocs(
     repo = "github.com/timmyfaraday/MultiStateSystems.jl.git"
)
