using Documenter
using GlobtimPlots

# docs/src/tutorial.md is derived from the README (see docs/tutorial_from_readme.jl);
# regenerate it on every build so the site can never drift from the README.
include("tutorial_from_readme.jl")
write_tutorial_page()

makedocs(
    sitename = "GlobtimPlots.jl",
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://gescholt.github.io/GlobtimPlots.jl",
        size_threshold = nothing,        # api.md (all docstrings) and the figure-bearing tutorial
        size_threshold_warn = nothing,
    ),
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
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
