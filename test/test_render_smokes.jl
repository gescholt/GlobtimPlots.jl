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
