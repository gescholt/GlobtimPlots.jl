# API Reference

```@meta
CurrentModule = GlobtimPlots
```

GlobtimPlots is organized by what you want to draw:

| Area | Functions |
|------|-----------|
| Level sets | [`cairo_plot_polyapprox_levelset`](@ref), [`plot_polyapprox_3d`](@ref) (GLMakie / WGLMakie), [`create_level_set_visualization`](@ref), [`create_level_set_animation`](@ref), [`plot_level_set`](@ref) |
| Morse diagnostics | [`plot_critical_eigenvalues`](@ref), [`plot_condition_numbers`](@ref), [`plot_hessian_norms`](@ref), [`plot_all_eigenvalues`](@ref), [`plot_raw_vs_refined_eigenvalues`](@ref) |
| Degree sweeps | [`plot_discrete_l2`](@ref), [`plot_convergence_analysis`](@ref), [`plot_l2_convergence_trials`](@ref), [`plot_distance_statistics`](@ref), [`plot_convergence_captured`](@ref), [`plot_theoretical_minimizer_distances`](@ref) |
| Subdivision | [`plot_subdivision_partition`](@ref), [`plot_subdivision_on_levelset`](@ref), [`plot_subdivision_on_levelset_from_bounds`](@ref), [`plot_subdivision_error_heatmap`](@ref), [`plot_subdivision_tree`](@ref), [`print_tree_summary`](@ref), [`plot_2d_partition`](@ref), [`plot_l2_trajectories`](@ref), [`interactive_error_explorer`](@ref), [`interactive_levelset_explorer`](@ref) |
| 1D approximation | [`plot_1d_polynomial_approximation`](@ref), [`plot_1d_comparison`](@ref) |
| Experiments & campaigns | [`create_experiment_plots`](@ref), [`create_campaign_comparison_plot`](@ref), [`create_single_plot`](@ref), [`generate_experiment_labels`](@ref), [`save_plot`](@ref), [`plot_experiment_results_static`](@ref), [`plot_experiment_results_interactive`](@ref) |
| Capture analysis | [`plot_capture_convergence`](@ref), [`plot_capture_sparsification_combined`](@ref), [`plot_refinement_trajectories!`](@ref) |
| Benchmark outcomes | [`plot_verdict_shares`](@ref), [`plot_predicate_pareto`](@ref) |

## Backends

GlobtimPlots draws through the abstract [Makie](https://docs.makie.org/) API, so every
function works with whichever backend is active. `CairoMakie` is a hard dependency, is
re-exported, and is active after `using GlobtimPlots`. The interactive backends are
optional and enable package extensions:

| Backend | Use for | Extension |
|---------|---------|-----------|
| `CairoMakie` | Static, publication-quality PNG / PDF / SVG | — (default) |
| `GLMakie` | Interactive windows, `plot_polyapprox_3d` with a recorded rotation, the explorers | `GlobtimPlotsGLMakieExt` |
| `WGLMakie` | Browser-embedded / Documenter-rendered interactive figures | `GlobtimPlotsWGLMakieExt` |

```julia
using GlobtimPlots           # CairoMakie active
fig = cairo_plot_polyapprox_levelset(apol, ainp, df_cp, df_min)
save("levelset.png", fig)

using GLMakie                # switches to GLMakie and loads the extension
plot_polyapprox_3d(apol, ainp, df_cp, df_min)
```

## Data interface

The level-set functions take the approximant and problem description through two small
abstract interfaces, [`AbstractPolynomialData`](@ref) and [`AbstractProblemInput`](@ref).
[`adapt_polynomial_data`](@ref) and [`adapt_problem_input`](@ref) wrap Globtim's
`ApproxPoly` and `TestInput`; [`GenericPolynomialData`](@ref) and
[`GenericProblemInput`](@ref) let you plot an approximation produced elsewhere. All
other functions are duck-typed over DataFrames and the result types of
GlobtimPostProcessing.

## Docstrings

```@autodocs
Modules = [GlobtimPlots]
```
