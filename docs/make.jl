using Documenter
using MultiStateSystems

makedocs(
    modules     = [MultiStateSystems],
    format      = Documenter.HTML(
        mathengine = Documenter.MathJax(),
        prettyurls = get(ENV, "CI", "false") == "true"
    ),
    sitename    = "MultiStateSystems.jl",
    authors     = "Tom Van Acker, Glenn Emmers",
    checkdocs   = :exports,  # Check that all exports are documented
    warnonly    = [:missing_docs, :cross_references],  # Only warn for these issues
    pages       = [ "Home"              => "index.md",
                    "Getting Started"   => "quickguide.md",
                    "DSL Manual"        =>
                        [ "State-Transition Diagram"=> "std.md",
                          "Distributions"           => "distribution.md",
                          "Generating Functions"    => "ugf.md",
                          "Network"                 => "network.md"],
                    "Models Manual"     =>
                        [ "Dependence"              => "dependence.md",
                          "Stochastic Processes"    => "processes.md",
                          "Generating Operators"    => "ugo.md"],
                    "Integrations"      =>
                        [ "DCIDE Integration"       => "dcide.md"],
                    "Output"            =>
                        [ "Indices"                 => "indices.md"],
                    "API Reference"     => "api_functions.md"
                  ]
)

deploydocs(
     repo = "github.com/timmyfaraday/MultiStateSystems.jl.git"
)
