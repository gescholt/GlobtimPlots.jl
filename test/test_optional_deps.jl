# Optional dependencies: GlobtimPlots needs no Makie backend, no GraphMakie and no
# GlobtimPostProcessing to load. Functions that need one of them are package
# extensions behind erroring stubs (src/extension_stubs.jl).
#
# The "not loaded" half runs in a fresh subprocess, because this test process has
# already loaded CairoMakie, GraphMakie and GlobtimPostProcessing.

using Test
using GlobtimPlots

function _run_subprocess(code::AbstractString)
    cmd = `$(Base.julia_cmd()) --startup-file=no --project=$(Base.active_project()) -e $code`
    return read(pipeline(cmd; stderr = devnull), String)
end

"Parse `KEY=value` lines printed by the subprocess."
_parse_kv(out) = Dict(
    String(m.captures[1]) => String(m.captures[2]) for
    m in eachmatch(r"^(\w+)=(.*)$"m, out)
)

@testset "Optional dependencies" begin
    @testset "bare `using GlobtimPlots`: no backend, no extensions" begin
        code = raw"""
        using GlobtimPlots
        ext(name) = Base.get_extension(GlobtimPlots, name) !== nothing
        errmsg(f) = try
            f(); "no error"
        catch e
            replace(sprint(showerror, e), '\n' => ' ')
        end
        println("CAIRO=", Base.find_package("CairoMakie") === nothing ? "absent" :
                          (Base.root_module_exists(Base.identify_package("CairoMakie")) ? "loaded" : "not_loaded"))
        println("EXT_GL=", ext(:GlobtimPlotsGLMakieExt))
        println("EXT_WGL=", ext(:GlobtimPlotsWGLMakieExt))
        println("EXT_GRAPH=", ext(:GlobtimPlotsGraphMakieExt))
        println("EXT_PP=", ext(:GlobtimPlotsPostProcessingExt))
        println("TREE=", errmsg(() -> plot_subdivision_tree(nothing)))
        println("CAPTURE=", errmsg(() -> plot_capture_convergence([], nothing)))
        println("SPARSE=", errmsg(() -> plot_capture_sparsification_combined([], nothing)))
        println("TRAJ=", errmsg(() -> plot_refinement_trajectories!(nothing, [], [])))
        fig = Figure(); lines(fig[1, 1], 1:3)
        println("SAVE=", errmsg(() -> save(tempname() * ".png", fig)))
        println("STYLE=", RefinementTrajectoryStyle() isa RefinementTrajectoryStyle)
        """
        kv = _parse_kv(_run_subprocess(code))
        @test get(kv, "CAIRO", "") == "not_loaded"
        @test get(kv, "EXT_GL", "") == "false"
        @test get(kv, "EXT_WGL", "") == "false"
        @test get(kv, "EXT_GRAPH", "") == "false"
        @test get(kv, "EXT_PP", "") == "false"
        @test occursin("requires GraphMakie", get(kv, "TREE", ""))
        for key in ("CAPTURE", "SPARSE", "TRAJ")
            @test occursin("requires GlobtimPostProcessing", get(kv, key, ""))
        end
        # Makie's own error when no backend is loaded
        @test occursin("No backend available", get(kv, "SAVE", ""))
        # Types stay in the core package
        @test get(kv, "STYLE", "") == "true"
    end

    @testset "`using GraphMakie` alone loads the GraphMakie extension" begin
        # Graphs and NetworkLayout are triggers too; they load with GraphMakie
        code = raw"""
        using GlobtimPlots, GraphMakie
        println("EXT_GRAPH=", Base.get_extension(GlobtimPlots, :GlobtimPlotsGraphMakieExt) !== nothing)
        println("EXT_PP=", Base.get_extension(GlobtimPlots, :GlobtimPlotsPostProcessingExt) !== nothing)
        """
        kv = _parse_kv(_run_subprocess(code))
        @test get(kv, "EXT_GRAPH", "") == "true"
        @test get(kv, "EXT_PP", "") == "false"
    end

    @testset "with the triggers loaded the features work" begin
        @eval using GraphMakie, GlobtimPostProcessing, CairoMakie
        @test Base.get_extension(GlobtimPlots, :GlobtimPlotsGraphMakieExt) !== nothing
        @test Base.get_extension(GlobtimPlots, :GlobtimPlotsPostProcessingExt) !== nothing

        tree = Base.invokelatest() do
            GlobtimPlots.Globtim.adaptive_refine(
                p -> sin(5 * p[1]) + 0.3 * sin(p[2]),
                [(-2.0, 2.0), (-2.0, 2.0)],
                4;
                l2_tolerance = 0.05,
                tolerance_mode = :absolute,
                max_depth = 3,
                parallel = false,
            )
        end
        fig = Base.invokelatest(plot_subdivision_tree, tree)
        @test fig isa Figure
        fig_std = Base.invokelatest(plot_subdivision_tree, tree; balanced_layout = false)
        @test fig_std isa Figure

        # Once the extension is loaded, a call matching none of its methods is an
        # ordinary MethodError, not the "requires" message
        @test_throws MethodError Base.invokelatest(plot_capture_convergence, 1, 2)
        @test_throws MethodError Base.invokelatest(plot_subdivision_tree)
    end
end
