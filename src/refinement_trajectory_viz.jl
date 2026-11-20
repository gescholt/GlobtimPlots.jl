"""
    Refinement Trajectory Visualization

Overlay local-refinement paths (Newton, trust-region, NelderMead) on an existing
Makie `Axis`.  Each trajectory is drawn from the raw polynomial critical point
(circle marker) to the refined result (triangle marker), with optional
intermediate iterates when the solver recorded a full trace.

The bang method `plot_refinement_trajectories!` mutates an existing axis so you
can compose it with any landscape / contour plot.

The style struct lives here; `plot_refinement_trajectories!` takes
GlobtimPostProcessing's `CriticalPointRefinementResult`, so its method is in the
GlobtimPostProcessing extension (`ext/GlobtimPlotsPostProcessingExt/`).
"""

# ── Style configuration ──────────────────────────────────────────────────────

"""
    RefinementTrajectoryStyle

Visual style for refinement trajectory overlays.

# Fields
- `linewidth::Float64`: Path line width (default 1.5)
- `markersize_start::Float64`: Start marker (circle) size (default 6.0)
- `markersize_end::Float64`: End marker (triangle) size (default 8.0)
- `marker_start::Symbol`: Marker for raw CP (default `:circle`)
- `marker_end::Symbol`: Marker for refined CP (default `:dtriangle`)
- `stroke_color`: Marker stroke color (default `:black`)
- `stroke_width::Float64`: Marker stroke width (default 0.5)
- `line_alpha::Float64`: Line opacity (default 0.8)
- `colors::Vector{Symbol}`: Palette for trajectory coloring (cycled)
"""
Base.@kwdef struct RefinementTrajectoryStyle
    linewidth::Float64 = 1.5
    markersize_start::Float64 = 6.0
    markersize_end::Float64 = 8.0
    marker_start::Symbol = :circle
    marker_end::Symbol = :dtriangle
    stroke_color::Any = :black
    stroke_width::Float64 = 0.5
    line_alpha::Float64 = 0.8
    colors::Vector{Symbol} = [
        :orangered,
        :darkorange,
        :gold,
        :yellowgreen,
        :cyan,
        :deepskyblue,
        :dodgerblue,
        :mediumpurple,
        :orchid,
        :hotpink,
    ]
end
