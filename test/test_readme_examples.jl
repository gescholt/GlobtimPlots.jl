# Executes the README tutorial.
#
# Every ```julia fenced block of README.md is extracted, written to a sandbox
# directory as readme_block_NN.jl, and included there in order inside one fresh
# module — the README states that later blocks reuse objects from earlier ones.
# Files the blocks save (`save("x.png", fig)`, `save_plot(fig, "x.png")`) must
# exist afterwards and be non-empty. ```julia-repl blocks are documentation
# only (GLMakie / Pkg-mode) and are not executed.

using Test
using DataFrames: nrow

const README_PATH = normpath(joinpath(@__DIR__, "..", "README.md"))

"Return the bodies of all ```julia fenced code blocks in `text`, in order."
function readme_julia_blocks(text::AbstractString)
    return [String(m.captures[1]) for m in eachmatch(r"```julia\n(.*?)```"s, text)]
end

"File names appearing in `save(\"...\")` / `save_plot(..., \"...\")` calls."
function readme_saved_files(text::AbstractString)
    return unique(
        String(m.captures[1]) for
        m in eachmatch(r"save\w*\([^\n]*?\"([^\"]+\.(?:png|pdf|svg))\"", text)
    )
end

@testset "README tutorial runs as written" begin
    readme = read(README_PATH, String)
    blocks = readme_julia_blocks(readme)
    @test length(blocks) >= 5

    saved = readme_saved_files(readme)
    @test !isempty(saved)

    sandbox = mktempdir()
    tutorial = Module(:ReadmeTutorial)
    cd(sandbox) do
        for (i, code) in enumerate(blocks)
            file = joinpath(sandbox, "readme_block_$(lpad(i, 2, '0')).jl")
            write(file, code)
            @testset "block $i" begin
                Base.include(tutorial, file)
                @test true
            end
        end
    end

    @testset "saved figure files" begin
        for name in saved
            path = joinpath(sandbox, name)
            @test isfile(path) && filesize(path) > 0
        end
    end

    # A few spot checks on the objects the tutorial builds
    @test isdefined(tutorial, :df_cp) && isdefined(tutorial, :results)
    @test nrow(tutorial.df_cp) > 0
    @test tutorial.labels == ["GN=100", "GN=200"]
end
