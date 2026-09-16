# The docs Tutorial page is derived from the README so there is one source of
# truth: test/test_readme_examples.jl executes the README's ```julia blocks, and
# the same blocks become Documenter @example blocks whose figures render at
# build time. docs/make.jl regenerates docs/src/tutorial.md on every build; the
# generated file is tracked so the mirror's file view has it, and the test suite
# checks that the tracked copy is in sync with the README.

"""
    tutorial_page_from_readme(readme_path) -> String

The Markdown of the docs Tutorial page: the README's `## Tutorial` section (up to
`## Function index`) with every ```julia fence turned into an `@example tutorial`
block. ```julia-repl fences (GLMakie, Pkg mode) are left as documentation.
"""
function tutorial_page_from_readme(readme_path::AbstractString)
    readme = read(readme_path, String)
    a = findfirst("## Tutorial\n", readme)
    b = findfirst("## Function index", readme)
    (a === nothing || b === nothing) && error("$readme_path: Tutorial section not found")
    section = readme[(last(a)+1):(first(b)-1)]
    section = replace(section, "```julia\n" => "```@example tutorial\n")
    return "# Tutorial\n\n*Generated from the package README by `docs/make.jl`; the test " *
           "suite executes the same code.*\n\n" *
           section
end

"Write the Tutorial page next to the other docs sources. Returns the path."
function write_tutorial_page(docs_dir::AbstractString = @__DIR__)
    path = joinpath(docs_dir, "src", "tutorial.md")
    write(path, tutorial_page_from_readme(joinpath(docs_dir, "..", "README.md")))
    return path
end
