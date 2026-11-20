# Stubs for functions whose methods live in package extensions.
#
# Each function has a catch-all method that raises an informative error when the
# extension is not loaded; the extension adds more specific methods, which
# dispatch selects, so nothing is overwritten. Once the extension is loaded, a
# call that matches none of its methods still lands here and gets the ordinary
# MethodError. The error text keeps the form "requires X.jl", like Globtim's
# HomotopyContinuation stubs.

function _require_extension(f, args, ext::Symbol, pkg::String)
    if Base.get_extension(@__MODULE__, ext) === nothing
        error(
            "$(nameof(f)) requires $(pkg).jl. Add `using $(pkg)` " *
            "(with $(pkg) in the active environment) to enable it.",
        )
    end
    throw(MethodError(f, args))
end

# --- GraphMakie extension (ext/GlobtimPlotsGraphMakieExt.jl) -------------------

"""
    plot_subdivision_tree(tree; kwargs...) -> Figure

Needs `using GraphMakie`: the method is in the GraphMakie package extension (its
triggers Graphs and NetworkLayout are GraphMakie dependencies and load with it).
Without it the call raises an error saying so.

Create publication-quality visualization of a subdivision tree.

# Arguments
- `tree`: SubdivisionTree from Globtim (uses duck typing)

# Keyword Arguments
- `style::TreeVizStyle = TreeVizStyle()`: Visual styling configuration
- `title::String = ""`: Optional figure title
- `show_legend::Bool = true`: Whether to show dimension color legend
- `balanced_layout::Bool = true`: Use split-position-aware layout (recommended)

# Returns
- `Figure`: Makie Figure; save it as PDF/SVG with CairoMakie loaded

# Notes
- Auto-scaling is enabled by default: node sizes, fonts, and figure dimensions
  adapt to tree size for optimal readability.
- Set `style=TreeVizStyle(auto_scale=false, auto_fig_size=false)` for fixed sizing.

# Example
```julia
using Globtim, GlobtimPlots, CairoMakie, GraphMakie

tree = adaptive_refine(f, bounds, 4)  # uses default relative L2 tolerance of 0.03
fig = plot_subdivision_tree(tree)
save("tree.pdf", fig)
```
"""
plot_subdivision_tree(args...; kwargs...) = _require_extension(
    plot_subdivision_tree,
    args,
    :GlobtimPlotsGraphMakieExt,
    "GraphMakie",
)

# --- GlobtimPostProcessing extension (ext/GlobtimPlotsPostProcessingExt/) ------

"""
    plot_capture_convergence(
        degree_capture_results::Vector{Tuple{Int, CaptureResult}},
        known::KnownCriticalPoints;
        fig_size = (900, 900),
        l2_errors = nothing,
        show_l2 = false,
        support_sizes = nothing,
        save_path = nothing,
    ) -> Figure

Needs `using GlobtimPostProcessing`: the method is in the GlobtimPostProcessing
package extension. Without it the call raises an error saying so.

Create a multi-panel figure showing capture analysis convergence.

When `support_sizes` is provided, the x-axis shows the number of nonzero polynomial
coefficients (support size) instead of polynomial degree. All spacing and bar widths
adapt automatically to the data via `_adaptive_spacing`.

# Panels (top to bottom)
1. **Distance distribution**: One sub-panel per CP type (min/max/saddle), each with its
   own y-range and log y-scale. Boxplots + jittered scatter of distances from each known
   CP to its nearest computed CP. Horizontal dashed lines at tolerance levels.
2. **Capture rate**: Grouped bars showing fraction of known CPs captured
   at each tolerance level.
3. **L2 approximation error** (optional): Only shown when both `l2_errors` data is
   provided AND `show_l2=true`. Not displayed by default.

# Arguments
- `degree_capture_results`: Vector of `(degree, CaptureResult)` tuples, one per degree.
  All CaptureResults must use the same `tolerance_fractions`.
- `known`: `KnownCriticalPoints` with the ground truth.

# Keyword Arguments
- `fig_size`: Figure dimensions in pixels.
- `l2_errors`: Optional `Vector{Tuple{Int, Float64}}` of `(degree, l2_error)` pairs.
- `show_l2`: Whether to display the L2 error panel. Default `false`. Requires `l2_errors`
  to be provided as well.
- `support_sizes`: Optional `Dict{Int, Int}` mapping degree -> support size (nonzero
  coefficients). When provided, x-axis shows support size instead of degree.
- `save_path`: Optional file path to save the figure (e.g., `"capture_convergence.png"`).

# Returns
- `Figure`: Makie figure object.

# Example
```julia
using GlobtimPlots, GlobtimPostProcessing, CairoMakie

# With support size on x-axis
support_sizes = Dict(dr.degree => dr.support_size
    for dr in degree_results if dr.status == "success")
fig = plot_capture_convergence(
    degree_capture_results, known_cps;
    support_sizes = support_sizes,
    l2_errors = [(dr.degree, dr.l2_approx_error) for dr in degree_results if dr.status == "success"],
)
display(fig)
```
"""
plot_capture_convergence(args...; kwargs...) = _require_extension(
    plot_capture_convergence,
    args,
    :GlobtimPlotsPostProcessingExt,
    "GlobtimPostProcessing",
)

"""
    plot_capture_sparsification_combined(
        entries,
        known::KnownCriticalPoints;
        fig_size = (1400, 700),
        save_path = nothing,
    ) -> Figure

Needs `using GlobtimPostProcessing`: the method is in the GlobtimPostProcessing
package extension. Without it the call raises an error saying so.

Two-panel side-by-side figure showing the sparsification speed vs accuracy tradeoff.

The x-axis uses support size (number of nonzero coefficients) from the Full variant
at each degree, so the axis reflects actual polynomial complexity. All spacing and
bar widths adapt automatically via `_adaptive_spacing`.

# Panel A (left) -- CP Accuracy: Distance to Known CPs
Stacked per-type subpanels with dodged boxplots for each support size x variant.
Log-y scale, tolerance reference lines. Shows how much CP accuracy degrades with
sparsification.

# Panel B (right) -- HC Solve Time
Grouped bar chart: x-axis = support size, bars dodged by variant. Y-axis = solve time (seconds).
Each sparsified bar is annotated with absolute time saved.
Shows how much faster the HC solve is with sparsification.

# Arguments
- `entries`: Vector of NamedTuples, each with fields:
  - `degree::Int` -- polynomial degree
  - `variant_label::String` -- "Full", "1e-5 (mild)", etc.
  - `threshold::Float64` -- sparsification threshold (0.0 for full polynomial)
  - `capture_result::CaptureResult` -- full capture analysis
  - `n_nonzero_coeffs::Int` -- number of non-zero coefficients
  - `l2_ratio::Float64` -- ratio of sparsified L2 norm to original (1.0 for full)
  - `solve_time::Float64` -- HC solve wall-clock time in seconds
- `known::KnownCriticalPoints`: Ground truth critical points.

# Keyword Arguments
- `fig_size`: Figure dimensions in pixels.
- `save_path`: Optional file path to save the figure.

# Returns
- `Figure`: Makie figure object.
"""
plot_capture_sparsification_combined(args...; kwargs...) = _require_extension(
    plot_capture_sparsification_combined,
    args,
    :GlobtimPlotsPostProcessingExt,
    "GlobtimPostProcessing",
)

"""
    plot_refinement_trajectories!(ax, refine_results, raw_points; kwargs...) -> Int

Needs `using GlobtimPostProcessing`: the method is in the GlobtimPostProcessing
package extension. Without it the call raises an error saying so.

Overlay refinement trajectories on an existing `Axis`.

Each trajectory visualizes how a local solver (Newton / trust-region / NelderMead)
moved a raw polynomial critical point toward the refined result.  Paths are drawn
as colored lines from circle (start) to triangle (end), with intermediate iterates
when the solver recorded a full trace.

# Arguments
- `ax::Axis`: Makie axis to draw on (mutated in place)
- `refine_results::Vector{CriticalPointRefinementResult}`: Per-CP refinement results
  (must have `.trace` field — a `Vector{Vector{Float64}}` of iterates)
- `raw_points::Vector{Vector{Float64}}`: Starting points (same length & ordering as
  `refine_results`)

# Keyword Arguments

## Selection — which trajectories to show
- `p_target::Union{Nothing, Vector{Float64}} = nothing`: Sort raw CPs by distance
  to this point and show the `n_show` closest.  Mutually exclusive with `indices`.
- `n_show::Int = 10`: Number of trajectories to show when using `p_target` sorting.
- `indices::Union{Nothing, Vector{Int}} = nothing`: Explicit 1-based indices into
  `refine_results` / `raw_points` to plot.  Overrides `p_target` / `n_show`.

## Projection — for N-D problems
- `dim_x::Int = 1`: Parameter dimension to map to the x-axis.
- `dim_y::Int = 2`: Parameter dimension to map to the y-axis.

## Clipping
- `clip_bounds::Union{Nothing, Vector{Vector{Float64}}} = nothing`: If provided,
  only draw trajectories that have at least one iterate inside this bounding box
  (format `[[lo1, hi1], [lo2, hi2], ...]`  — checked on `dim_x` and `dim_y` only).

## Visual style
- `style::RefinementTrajectoryStyle = RefinementTrajectoryStyle()`: Full style config.
- `label::Union{String, Nothing} = "refinement paths"`: Legend label for the first
  trajectory group.  Pass `nothing` to suppress legend entries.

# Returns
- `Int`: Number of trajectories actually drawn.

# Example
```julia
using GlobtimPlots, GlobtimPostProcessing, CairoMakie

fig = Figure()
ax = Axis(fig[1,1])
contourf!(ax, x_range, y_range, Z; colormap=:inferno)

n_drawn = plot_refinement_trajectories!(ax, results, raw_cps;
    p_target = p_true, n_show = 15,
    style = RefinementTrajectoryStyle(linewidth=2.0))
```
"""
plot_refinement_trajectories!(args...; kwargs...) = _require_extension(
    plot_refinement_trajectories!,
    args,
    :GlobtimPlotsPostProcessingExt,
    "GlobtimPostProcessing",
)
