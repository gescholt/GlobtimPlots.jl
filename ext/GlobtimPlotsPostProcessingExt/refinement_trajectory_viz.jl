
# Refinement trajectory overlay (GlobtimPostProcessing extension). The style
# struct `RefinementTrajectoryStyle` is in src/refinement_trajectory_viz.jl and the
# docstring on the stub in src/extension_stubs.jl.

function GlobtimPlots.plot_refinement_trajectories!(
    ax,
    refine_results::Vector{CriticalPointRefinementResult},
    raw_points::Vector{Vector{Float64}};
    # Selection
    p_target::Union{Nothing,Vector{Float64}} = nothing,
    n_show::Int = 10,
    indices::Union{Nothing,Vector{Int}} = nothing,
    # Projection
    dim_x::Int = 1,
    dim_y::Int = 2,
    # Clipping
    clip_bounds::Union{Nothing,Vector{Vector{Float64}}} = nothing,
    # Style
    style::RefinementTrajectoryStyle = RefinementTrajectoryStyle(),
    label::Union{String,Nothing} = "refinement paths",
)::Int
    length(refine_results) == length(raw_points) || error(
        "refine_results ($(length(refine_results))) and raw_points ($(length(raw_points))) must have the same length",
    )

    # ── Determine which trajectories to plot ─────────────────────────────────
    plot_indices = if indices !== nothing
        # Explicit selection
        for idx in indices
            1 <= idx <= length(refine_results) ||
                error("index $idx out of range [1, $(length(refine_results))]")
        end
        indices
    elseif p_target !== nothing
        # Sort by distance to target, take top n_show
        dists = [norm(rp .- p_target) for rp in raw_points]
        order = sortperm(dists)
        order[1:min(n_show, length(order))]
    else
        # No selection criteria — show first n_show
        collect(1:min(n_show, length(refine_results)))
    end

    # ── Draw trajectories ────────────────────────────────────────────────────
    n_drawn = 0
    legend_added = false
    colors = style.colors

    for (rank, idx) in enumerate(plot_indices)
        rr = refine_results[idx]
        tpts = rr.trace
        length(tpts) < 2 && continue

        tx = [Float64(p[dim_x]) for p in tpts]
        ty = [Float64(p[dim_y]) for p in tpts]

        # Clip check: skip if no iterate falls inside the clip window
        if clip_bounds !== nothing
            xb = clip_bounds[1]  # [lo, hi] for dim_x
            yb = clip_bounds[2]  # [lo, hi] for dim_y
            any_inside = any(
                xb[1] <= tx[k] <= xb[2] && yb[1] <= ty[k] <= yb[2] for k in eachindex(tx)
            )
            any_inside || continue
        end

        col = colors[mod1(rank, length(colors))]

        # Legend entry only on first drawn trajectory
        lbl = (!legend_added && label !== nothing) ? label : nothing

        lines!(
            ax,
            tx,
            ty;
            color = (col, style.line_alpha),
            linewidth = style.linewidth,
            label = lbl,
        )

        # Start marker (raw CP)
        scatter!(
            ax,
            [tx[1]],
            [ty[1]];
            marker = style.marker_start,
            markersize = style.markersize_start,
            color = col,
            strokecolor = style.stroke_color,
            strokewidth = style.stroke_width,
        )

        # End marker (refined CP)
        scatter!(
            ax,
            [tx[end]],
            [ty[end]];
            marker = style.marker_end,
            markersize = style.markersize_end,
            color = col,
            strokecolor = style.stroke_color,
            strokewidth = style.stroke_width,
        )

        legend_added = true
        n_drawn += 1
    end

    return n_drawn
end
