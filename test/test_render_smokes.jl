# Render smokes for the most-used GlobtimPlots entry points.
#
# Coverage targets picked by repo-wide call-site count:
#   cairo_plot_polyapprox_levelset (10 call sites)
#   create_level_set_visualization  (5 call sites)
#   create_experiment_plots         (4 call sites)
# plot_polyapprox_3d is GLMakie-extension-only (needs a GL context) and is
# exercised by the interactive explorers, not here.

using Test
using GlobtimPlots
using CairoMakie
using DataFrames
using StaticArrays

@testset "Render smokes for top entry points" begin
    @testset "cairo_plot_polyapprox_levelset" begin
        # Synthetic 2D tensor grid with z = x1^2 + x2^2
        n = 21
        xs = range(-1.0, 1.0; length = n)
        grid = Matrix{Float64}(undef, n * n, 2)
        z = Vector{Float64}(undef, n * n)
        for (k, (x, y)) in enumerate(Iterators.product(xs, xs))
            grid[k, 1] = x
            grid[k, 2] = y
            z[k] = x^2 + y^2
        end
        pol = GenericPolynomialData(coeffs = nothing, grid = grid, z = z)
        TR = GenericProblemInput(dim = 2, center = [0.0, 0.0], sample_range = 1.0)

        df = DataFrame(
            x1 = [0.5, -0.5, 0.2],
            x2 = [0.5, -0.5, -0.3],
            z = [0.5, 0.5, 0.13],
            close = [true, false, true],
        )
        df_min = DataFrame(
            x1 = [0.0, 0.7],
            x2 = [0.0, 0.7],
            value = [0.0, 0.98],
            captured = [true, false],
        )

        fig = cairo_plot_polyapprox_levelset(pol, TR, df, df_min)
        @test fig isa Figure

        # Renders to a non-empty PNG
        tmpfile = tempname() * ".png"
        save(tmpfile, fig)
        @test isfile(tmpfile) && filesize(tmpfile) > 0
        rm(tmpfile; force = true)

        # Chebyshev level spacing path
        fig2 = cairo_plot_polyapprox_levelset(pol, TR, df, df_min; chebyshev_levels = true)
        @test fig2 isa Figure
    end

    @testset "create_level_set_visualization (3D grid)" begin
        n = 9
        xs = range(-1.0, 1.0; length = n)
        grid = [SVector{3,Float64}(x, y, z) for x in xs, y in xs, z in xs]
        f = p -> p[1]^2 + p[2]^2 + p[3]^2

        result = create_level_set_visualization(f, grid, nothing, (0.2, 0.8))
        @test result !== nothing
    end

    @testset "level sets from values, with index-coloured overlays" begin
        n = 7
        xs = range(-1.0, 1.0; length = n)
        grid = [SVector{3,Float64}(x, y, z) for x in xs, y in xs, z in xs]
        vals = [p[1]^2 - p[2]^2 + p[3]^2 for p in grid]
        truth = DataFrame(x1 = [0.0], x2 = [0.0], x3 = [0.0], z = [0.0], index = [1])
        cands = DataFrame(x1 = [0.01, 0.5], x2 = [0.0, 0.5], x3 = [0.0, 0.0], z = [0.0, 0.0], index = [1, 0])
        ov = [(; df = truth, label = "truth", marker = :circle), (; df = cands, label = "surrogate", marker = :diamond)]
        fig = create_level_set_visualization(vals, grid, nothing, (-0.5, 0.5); overlays = ov, title = "saddle")
        @test fig isa Figure
        mp4 = tempname() * ".mp4"
        create_level_set_animation(vals, grid, nothing, (-0.5, 0.5); overlays = ov, fps = 2, duration = 1, output_file = mp4)
        @test isfile(mp4) && filesize(mp4) > 0
        rm(mp4; force = true)
        # edge crossings: on the level set to within one cell's variation of f, drawn the same way
        xs_c = level_set_crossings(grid, vals, 0.25)
        h = step(xs)
        @test !isempty(xs_c)
        @test all(p -> abs(p[1]^2 - p[2]^2 + p[3]^2 - 0.25) < 3h * (2 + h), xs_c)
        @test isempty(level_set_crossings(grid, vals, 10.0))
        create_level_set_animation(vals, grid, nothing, (-0.5, 0.5); overlays = ov, fps = 2, duration = 1,
            output_file = mp4, level_style = :crossings)
        @test isfile(mp4) && filesize(mp4) > 0
        rm(mp4; force = true)
        @test_throws ArgumentError create_level_set_animation(vals, grid, nothing, (-0.5, 0.5); level_style = :bogus)
        @test_throws ArgumentError create_level_set_visualization(vals, grid, nothing, (-0.5, 0.5);
            overlays = [(; df = DataFrame(x1 = [0.0]), label = "bad")])
    end

    @testset "create_sublevel_animation (two wells merging)" begin
        n = 11
        xs = range(-1.0, 1.0; length = n)
        grid = [SVector{3,Float64}(x, y, z) for x in xs, y in xs, z in xs]
        vals = [(p[1]^2 - 0.5)^2 + p[2]^2 + p[3]^2 for p in grid]
        saddle = DataFrame(x1 = [0.0], x2 = [0.0], x3 = [0.0], z = [0.25], index = [1])
        mp4 = tempname() * ".mp4"
        fig = create_sublevel_animation(vals, grid, (0.0, 0.6); fps = 2, duration = 2,
            output_file = mp4, overlays = [(; df = saddle, label = "saddle")])
        @test fig isa Figure
        @test isfile(mp4) && filesize(mp4) > 0
        rm(mp4; force = true)
        # two wells at x₁ = ±√0.5 (value 0) merging at the saddle (value 0.25 on the grid)
        bars = sublevel_barcode(vals)
        @test length(bars) == 2
        @test count(b -> isinf(b.death), bars) == 1
        young = only(filter(b -> isfinite(b.death), bars))
        @test 0.2 <= young.death <= 0.3
    end

    @testset "create_level_set_sweep_animation (a fold: two critical points appear)" begin
        n = 9
        xs = range(-1.0, 1.0; length = n)
        grid = [SVector{3,Float64}(x, y, z) for x in xs, y in xs, z in xs]
        ts = (-0.3, 0.0, 0.3)
        frames = [[p[1]^3 - t * p[1] + p[2]^2 + p[3]^2 for p in grid] for t in ts]
        cps(t) = t <= 0 ? DataFrame(x1 = Float64[], x2 = Float64[], x3 = Float64[], z = Float64[], index = Int[]) :
                 DataFrame(x1 = [-sqrt(t / 3), sqrt(t / 3)], x2 = [0.0, 0.0], x3 = [0.0, 0.0], z = [0.0, 0.0], index = [1, 0])
        ov = [[(; df = cps(t), label = "critical points")] for t in ts]
        mp4 = tempname() * ".mp4"
        fig = create_level_set_sweep_animation(frames, grid, 0.2,
            VisualizationParameters(point_tolerance = 0.05, point_window = 1e9);
            frame_overlays = ov, frame_labels = ["t = $t" for t in ts], substeps = 2, fps = 4, output_file = mp4)
        @test fig isa Figure
        @test isfile(mp4) && filesize(mp4) > 0
        rm(mp4; force = true)
        @test_throws ArgumentError create_level_set_sweep_animation(frames, grid, 0.2; frame_labels = ["a"])
    end

    @testset "create_experiment_plots (Static backend)" begin
        result =
            (experiment_id = "render_smoke", enabled_tracking = ["approximation_quality"])
        stats = Dict(
            "approximation_quality" =>
                Dict("degrees" => [4, 6, 8], "l2_errors" => [0.5, 0.1, 0.02]),
        )

        fig = create_experiment_plots(result, stats; backend = Static)
        @test fig isa Figure

        tmpfile = tempname() * ".png"
        save(tmpfile, fig)
        @test isfile(tmpfile) && filesize(tmpfile) > 0
        rm(tmpfile; force = true)
    end
end

@testset "Benchmark outcome plots (qjf0)" begin
    @testset "plot_verdict_shares" begin
        df = DataFrame(
            family = ["lv4d", "lv4d", "lv4d", "fhn3d", "fhn3d", "fhn3d"],
            winner = ["per_axis", "global", "tie", "per_axis", "per_axis", "global"],
            base_degree = [2, 2, 3, 2, 2, 2],
            budget = [4, 4, 8, 4, 8, 8],
        )
        fig = plot_verdict_shares(df)
        @test fig isa Figure
        fig2 = plot_verdict_shares(df; strat_cols = [:budget])
        @test fig2 isa Figure
        @test_throws ErrorException plot_verdict_shares(first(df, 0))
    end

    @testset "plot_predicate_pareto" begin
        df = DataFrame(
            trial = ["t1", "t2", "t3", "t4"],
            mean_evals = [100.0, 200.0, 150.0, 400.0],
            mean_recovery = [1e-2, 1e-4, 5e-3, 1e-4],
            pareto_score = [0.5, 0.9, 0.6, 0.4],
        )
        fig = plot_predicate_pareto(df; label_col = :trial)
        @test fig isa Figure
        # precomputed front column is honored verbatim
        df.is_pareto = [true, true, false, false]
        fig2 = plot_predicate_pareto(df; front_col = :is_pareto, score_col = nothing)
        @test fig2 isa Figure
        @test_throws ErrorException plot_predicate_pareto(first(df, 0))
    end
end
