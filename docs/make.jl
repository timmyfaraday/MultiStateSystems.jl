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
                          "Generating Functions"    => "ugf.md",
                          "Network"                 => "network.md"],
                    "Models Manual"     =>
                        [ "Stochastic Processes"    => "processes.md"
                          "Generating Operators"    => "ugo.md"]
                  ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
     repo = "github.com/timmyfaraday/MultiStateSystems.jl.git"
)
