using Documenter
using GlobtimPlots

makedocs(
    sitename = "GlobtimPlots.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://gescholt.github.io/GlobtimPlots.jl",
    ),
    pages = [
        "Home" => "index.md",
        "API Reference" => "api.md",
    ],
    modules = [GlobtimPlots],
)

deploydocs(
    repo = "github.com/gescholt/GlobtimPlots.jl.git",
    target = "build",
    branch = "gh-pages",
    devbranch = "main",
)
