module GlobtimPlots

# Core plotting dependencies
using CairoMakie
# NOTE: GLMakie is not a hard dependency - it causes a precompilation segfault on
# macOS (GLFW monitor detection bug). Load it separately: using GLMakie; GLMakie.activate!()
using DataFrames
using Statistics
using ColorSchemes

# Re-export CairoMakie so `using GlobtimPlots` alone gives a working static backend
# (Figure, Axis, save, ...). Users can still switch: using GLMakie; GLMakie.activate!()
using Reexport
@reexport using CairoMakie

# Data structures from Globtim / GlobtimPostProcessing are consumed by duck typing
# (plus the adapters in interfaces.jl), which keeps the dependency direction one-way.

# Abstract interfaces + adapters for Globtim's ApproxPoly / TestInput
include("interfaces.jl")

# Core plotting functionality
include("analysis_plots.jl")            # level sets, convergence, distance statistics
include("level_set_viz.jl")             # 3D level set visualization
include("experiment_results_plots.jl")  # cluster experiment results visualization

# Campaign plotting (duck-typed over ExperimentResult / CampaignResults)
include("CampaignPlotting.jl")

# Subdivision tree visualization (adaptive refinement from Globtim)
include("subdivision_tree_viz.jl")

# 1D polynomial approximation visualization
include("plot_1d_approx.jl")

# 2D partition visualization
# (GeometryBasics is a transitive dep via CairoMakie)
using GeometryBasics: Rect
include("partition_viz.jl")

# Subdivision partition visualization (adaptive subdivision experiments)
include("subdivision_partition_viz.jl")

# Subdivision error heatmap (runtime plot: requires tree with polynomials in memory)
include("subdivision_error_heatmap.jl")

# Interactive explorers (zoom-refine, click-to-optimize; need GLMakie for interactivity)
include("interactive_error_explorer.jl")
include("interactive_levelset_explorer.jl")

# Capture analysis + refinement trajectories (GlobtimPostProcessing types)
include("capture_analysis_plots.jl")
include("refinement_trajectory_viz.jl")

# Backend-agnostic Hessian/eigenvalue plots
include("hessian_eigenvalue_plots.jl")

# Backend-agnostic 3D polynomial approximation surface (GLMakie extension adds record())
include("polyapprox_3d.jl")

# Benchmark-outcome plots: counterfactual verdict shares + predicate Pareto
include("benchmark_outcome_plots.jl")

# Abstract types and adapters
export AbstractPolynomialData, AbstractProblemInput, AbstractCriticalPointData
export GenericPolynomialData, GenericProblemInput
export adapt_polynomial_data, adapt_problem_input
export _plot_polyapprox_3d_impl

# Level sets, convergence and distance statistics (analysis_plots.jl)
export cairo_plot_polyapprox_levelset, plot_convergence_analysis, plot_discrete_l2
export plot_distance_statistics, plot_convergence_captured
export plot_theoretical_minimizer_distances
export plot_l2_convergence_trials
export plot_experiment_results_static, plot_experiment_results_interactive

"""
    plot_polyapprox_3d(pol, TR, df, df_min; kwargs...) -> Figure

3D surface of the polynomial approximant `pol` over the domain described by `TR`,
with the critical points of `df` and the refined minima of `df_min` drawn on it.
`pol` and `TR` implement [`AbstractPolynomialData`](@ref) / [`AbstractProblemInput`](@ref)
(see [`adapt_polynomial_data`](@ref), [`adapt_problem_input`](@ref)).

The methods live in the backend extensions: load `GLMakie` (adds `rotate = true` /
`filename` for a recorded turntable animation) or `WGLMakie` (WebGL widget) first.
Without one of them the call raises a `MethodError`. The backend-agnostic drawing
code is `_plot_polyapprox_3d_impl` in `src/polyapprox_3d.jl`.

# Keyword arguments
- `figure_size = (1000, 800)`, `z_limits = nothing`, `alpha_surface = 0.7`
- `show_captured = true`: draw the refined minima of `df_min`
- `fade = false`, `z_cut = 0.25`: fade the surface above a value cut
- `color_by = :proximity`: marker coloring of the critical points
"""
function plot_polyapprox_3d end
export plot_polyapprox_3d

# Hessian / eigenvalue diagnostics (hessian_eigenvalue_plots.jl, backend-agnostic)
export plot_hessian_norms, plot_condition_numbers, plot_critical_eigenvalues
export plot_all_eigenvalues, plot_raw_vs_refined_eigenvalues

# Campaign plotting (CampaignPlotting.jl)
export PlotBackend, Interactive, Static
export create_experiment_plots,
    create_campaign_comparison_plot, create_single_plot, save_plot
export generate_experiment_labels

# Level set visualization (level_set_viz.jl)
export LevelSetData, VisualizationParameters
export prepare_level_set_data, to_makie_format, plot_level_set
export create_level_set_visualization, create_level_set_animation

# Utility functions
export transform_coordinates, points_in_hypercube
export analyze_convergence_distances, analyze_captured_distances

# Subdivision tree visualization
export plot_subdivision_tree, TreeVizStyle, print_tree_summary

# 1D polynomial approximation visualization
export plot_1d_polynomial_approximation, plot_1d_comparison

# 2D partition visualization
export plot_2d_partition, plot_l2_trajectories, PartitionStyle

# Subdivision partition visualization
export plot_subdivision_partition,
    plot_subdivision_on_levelset,
    plot_subdivision_on_levelset_from_bounds,
    SubdivisionPartitionStyle
export format_degree_label, effective_degree, min_degree
export CP_TYPE_TABLE, cp_type_appearance, normalize_cp_key

# Subdivision error heatmap
export plot_subdivision_error_heatmap, find_containing_leaf, compute_error_grid

# Interactive explorers
export interactive_error_explorer
export interactive_levelset_explorer

# Capture analysis visualization
export plot_capture_convergence, plot_capture_sparsification_combined

# Refinement trajectory visualization
export plot_refinement_trajectories!, RefinementTrajectoryStyle

# Benchmark-outcome plots
export plot_verdict_shares, plot_predicate_pareto

end # module GlobtimPlots
