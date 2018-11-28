# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [fixdoctests]
#
# for local builds.

using Documenter
using LegendHDF5IO

makedocs(
    sitename = "LegendHDF5IO",
    modules = [LegendHDF5IO],
    format = :html,
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    html_prettyurls = !("local" in ARGS),
    html_canonical = "https://legend-exp.github.io/LegendHDF5IO.jl/stable/",
)

deploydocs(
    repo = "github.com/legend-exp/LegendHDF5IO.jl.git"
)
