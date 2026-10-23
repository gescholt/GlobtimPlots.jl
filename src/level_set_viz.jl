"""
Level Set Visualization Functions for GlobtimPlots

3D level set visualization (Axis3 scatter of near-level grid points, with a level
slider). Backend-agnostic Makie code: static under CairoMakie, interactive once
GLMakie is activated.

LevelSetData and VisualizationParameters types are defined in Globtim
and re-exported here for convenience.
"""

# NOTE: GLMakie removed from direct import - causes macOS precompilation segfault
# Users who need 3D interactive plots should: using GLMakie; GLMakie.activate!()
using Makie
using StaticArrays
using DataFrames

# Import canonical type definitions from Globtim
using Globtim: LevelSetData, VisualizationParameters

# ── point-set overlays (precomputed-values overloads) ─────────────────────────
#
# `overlays` is a vector of NamedTuples `(; df, label, marker = :circle,
# markersize = 16, color = :darkorange)`. `df` has columns x1, x2, x3 (the point)
# and z (its value on the scale of `values`); with an `index` column the points are
# coloured by Morse index (INDEX_COLORS, 0 = minimum) and `color` is ignored. A
# point is drawn while |z − level| ≤ `params.point_window`. One scatter per
# (overlay, index), each with a fixed colour, so a level change only replaces
# point lists.

const INDEX_COLORS = (:royalblue, :darkorange, :crimson, :purple, :black)
const INDEX_NAMES = ("minimum", "index 1", "index 2", "index 3", "index ≥ 4")

function _overlay_layers(ax, overlays)
    layers = []
    for ov in overlays
        df = ov.df
        all(c -> hasproperty(df, c), (:x1, :x2, :x3, :z)) ||
            throw(ArgumentError("overlay $(ov.label): df needs columns x1, x2, x3, z"))
        marker = get(ov, :marker, :circle)
        ms = get(ov, :markersize, 16)
        if hasproperty(df, :index)
            for k in sort(unique(min.(df.index, 4)))
                k >= 0 || continue
                obs = Observable(Point3f[])
                scatter!(
                    ax,
                    obs;
                    color = INDEX_COLORS[k+1],
                    marker = marker,
                    markersize = ms,
                    strokewidth = 1,
                    strokecolor = :black,
                )
                push!(layers, (; rows = df[min.(df.index, 4).==k, :], obs = obs))
            end
        else
            obs = Observable(Point3f[])
            scatter!(
                ax,
                obs;
                color = get(ov, :color, :darkorange),
                marker = marker,
                markersize = ms,
                strokewidth = 1,
                strokecolor = :black,
            )
            push!(layers, (; rows = df, obs = obs))
        end
    end
    return layers
end

function _update_overlays!(layers, level, window)
    for L in layers
        L.obs[] = Point3f[
            Point3f(r.x1, r.x2, r.x3) for r in eachrow(L.rows) if
            !any(isnan, (r.x1, r.x2, r.x3, r.z)) && abs(r.z - level) <= window
        ]
    end
end

# legend: one grey marker per overlay (its source), one colour per Morse index present
function _overlay_legend!(pos, overlays; index_names = INDEX_NAMES)
    elems = Any[]
    labels = String[]
    for ov in overlays
        push!(
            elems,
            MarkerElement(
                marker = get(ov, :marker, :circle),
                color = :grey70,
                strokecolor = :black,
                strokewidth = 1,
                markersize = 14,
            ),
        )
        push!(labels, ov.label)
    end
    idx = sort(
        unique(
            reduce(
                vcat,
                [min.(ov.df.index, 4) for ov in overlays if hasproperty(ov.df, :index)];
                init = Int[],
            ),
        ),
    )
    for k in idx
        k >= 0 || continue
        push!(
            elems,
            MarkerElement(marker = :circle, color = INDEX_COLORS[k+1], markersize = 14),
        )
        push!(labels, index_names[k+1])
    end
    Legend(pos, elems, labels; framevisible = false)
end

"""
    prepare_level_set_data(
        grid::Array{SVector{3,T},3},
        values::Array{T},
        level::T;
        tolerance::T = T(0.01)
    ) where {T<:AbstractFloat}

Identify points in the grid that are close to the specified level.

# Arguments
- `grid`: 3D array of grid points
- `values`: Function values at each grid point
- `level`: Target level value
- `tolerance`: Maximum distance from level to include point

# Returns
- `LevelSetData{T}`: Container with filtered points and values
"""
function prepare_level_set_data(
    grid::Array{SVector{3,T},3},
    values::Array{T},
    level::T;
    tolerance::T = T(0.01),
) where {T<:AbstractFloat}
    size(grid) == size(values) ||
        throw(DimensionMismatch("Grid and values must have same dimensions"))
    tolerance > zero(T) || throw(ArgumentError("Tolerance must be positive"))

    flat_grid = vec(grid)
    flat_values = vec(values)

    level_set_mask = @. abs(flat_values - level) < tolerance

    return LevelSetData(flat_grid[level_set_mask], flat_values[level_set_mask], level)
end

"""
    to_makie_format(level_set::LevelSetData{T}) where {T<:AbstractFloat}

Convert LevelSetData to format suitable for Makie plotting.

# Returns
Named tuple with:
- `points`: 3×N matrix of coordinates
- `values`: vector of function values
- `xyz`: tuple of (x, y, z) coordinate vectors
"""
function to_makie_format(level_set::LevelSetData{T}) where {T<:AbstractFloat}
    isempty(level_set.points) &&
        return (points = Matrix{T}(undef, 3, 0), values = T[], xyz = (T[], T[], T[]))

    points = reduce(hcat, level_set.points)
    return (
        points = points,
        values = level_set.values,
        xyz = (view(points, 1, :), view(points, 2, :), view(points, 3, :)),
    )
end

"""
    plot_level_set(
        formatted_data;
        fig_size = (800, 600),
        marker_size = 4,
        title = "Level Set Visualization"
    )

Create a basic 3D scatter plot of level set points.

# Arguments
- `formatted_data`: Output from `to_makie_format`
- `fig_size`: Figure dimensions
- `marker_size`: Size of scatter markers
- `title`: Plot title

# Returns
- GLMakie Figure object
"""
function plot_level_set(
    formatted_data;
    fig_size = (800, 600),
    marker_size = 4,
    title = "Level Set Visualization",
)
    fig = Figure(size = fig_size)
    ax = Axis3(fig[1, 1], title = title, xlabel = "x₁", ylabel = "x₂", zlabel = "x₃")

    scatter!(ax, formatted_data.xyz..., markersize = marker_size)

    display(fig)
    return fig
end

"""
    create_level_set_visualization(
        f,
        grid::Array{SVector{3,T},3},
        df::Union{DataFrame,Nothing},
        z_range::Tuple{T,T},
        params::VisualizationParameters{T} = VisualizationParameters{T}()
    ) where {T<:AbstractFloat}

Create an interactive 3D level set visualization with a slider.

# Arguments
- `f`: Objective function to visualize
- `grid`: 3D array of sample points
- `df`: Optional DataFrame with critical points (columns: x1, x2, x3, z)
- `z_range`: Tuple (min, max) for level values
- `params`: Visualization parameters

# Returns
- GLMakie Figure with interactive slider

# Example
```julia
using GLMakie, StaticArrays, GlobtimPlots

# Create a simple test function
f(x) = x[1]^2 + x[2]^2 + x[3]^2

# Create grid
grid = [SVector{3,Float64}(x, y, z)
        for x in range(-1, 1, 50),
            y in range(-1, 1, 50),
            z in range(-1, 1, 50)]

# Visualize
fig = create_level_set_visualization(f, grid, nothing, (0.0, 3.0))
```
"""
function create_level_set_visualization(
    f,
    grid::Array{SVector{3,T},3},
    df::Union{DataFrame,Nothing},
    z_range::Tuple{T,T},
    params::VisualizationParameters{T} = VisualizationParameters{T}(),
) where {T<:AbstractFloat}
    grid_points = vec(grid)
    valid_points = filter(p -> !any(isnan, p), grid_points)
    isempty(valid_points) && throw(ArgumentError("Grid contains no valid points"))

    z_min, z_max = z_range
    (isnan(z_min) || isnan(z_max)) && throw(ArgumentError("Invalid z_range"))

    fig = Figure(size = params.fig_size)
    ax = Axis3(fig[1, 1], xlabel = "x₁", ylabel = "x₂", zlabel = "x₃")

    x_range = extrema(p[1] for p in valid_points)
    y_range = extrema(p[2] for p in valid_points)
    z_range_grid = extrema(p[3] for p in valid_points)

    limits!(ax, x_range..., y_range..., z_range_grid...)

    level_slider =
        Slider(fig[2, 1], range = range(z_min, z_max, length = 1000), startvalue = z_min)

    Label(
        fig[3, 1],
        @lift(string("Level: ", round($(level_slider.value), digits = 3))),
        tellwidth = false,
    )

    level_points = Observable(Point3f[])
    data_points = Observable(Point3f[])

    # Pre-compute function values for the entire grid
    values = zeros(T, size(grid)...)
    @inbounds for i in eachindex(grid_points)
        point = grid_points[i]
        values[i] = any(isnan, point) ? NaN : f(point)
    end

    scatter!(
        ax,
        level_points,
        color = :blue,
        markersize = 6,
        alpha = 0.7,
        label = "Level Set",
    )

    if !isnothing(df)
        scatter!(
            ax,
            data_points,
            color = :darkorange,
            marker = :diamond,
            markersize = 30,
            label = "Data Points",
        )
    end

    function update_visualization(level::T) where {T<:AbstractFloat}
        try
            # Update level set points
            level_data = prepare_level_set_data(
                grid,
                values,
                level,
                tolerance = params.point_tolerance,
            )

            formatted_data = to_makie_format(level_data)

            # Update points atomically
            new_points = Point3f[]
            if !isempty(formatted_data.xyz[1])
                for (x, y, z) in zip(formatted_data.xyz...)
                    if !any(isnan, (x, y, z))
                        push!(new_points, Point3f(x, y, z))
                    end
                end
            end
            level_points[] = new_points

            if !isnothing(df)
                visible_points = Point3f[]
                for row in eachrow(df)
                    if !any(isnan, [row.x1, row.x2, row.x3, row.z]) &&
                       abs(row.z - level) ≤ params.point_tolerance
                        push!(visible_points, Point3f(row.x1, row.x2, row.x3))
                    end
                end
                data_points[] = visible_points
            end
        catch e
            @error "Error in visualization update" exception = e
            rethrow(e)
        end
    end

    on(level_slider.value) do level
        update_visualization(level)
    end

    update_visualization(z_min)
    axislegend(ax, position = :rt)

    return fig
end

"""
    create_level_set_animation(
        f,
        grid::Array{SVector{3,T},3},
        df::Union{DataFrame,Nothing},
        z_range::Tuple{T,T},
        params::VisualizationParameters{T} = VisualizationParameters{T}();
        fps::Int = 30,
        duration::Int = 20,
        output_file::String = "level_set_animation.mp4"
    ) where {T<:AbstractFloat}

Create an animated level set visualization that sweeps through levels
while rotating the camera.

# Arguments
- `f`: Objective function
- `grid`: 3D array of sample points
- `df`: Optional DataFrame with critical points
- `z_range`: Range of level values to animate
- `params`: Visualization parameters
- `fps`: Frames per second for animation
- `duration`: Animation duration in seconds
- `output_file`: Path to save the animation

# Returns
- GLMakie Figure object

# Example
```julia
fig = create_level_set_animation(
    f, grid, df, (-3.0, 6.0);
    fps = 30,
    duration = 10,
    output_file = "my_animation.mp4"
)
```
"""
function create_level_set_animation(
    f,
    grid::Array{SVector{3,T},3},
    df::Union{DataFrame,Nothing},
    z_range::Tuple{T,T},
    params::VisualizationParameters{T} = VisualizationParameters{T}();
    fps::Int = 30,
    duration::Int = 20,
    output_file::String = "level_set_animation.mp4",
) where {T<:AbstractFloat}
    grid_points = vec(grid)
    valid_points = filter(p -> !any(isnan, p), grid_points)

    z_min, z_max = z_range

    fig = Figure(size = params.fig_size)
    ax = Axis3(
        fig[1, 1],
        title = "Level Set Visualization",
        xlabel = "x₁",
        ylabel = "x₂",
        zlabel = "x₃",
    )

    # Set up initial ranges and limits
    x_range = extrema(p[1] for p in valid_points)
    y_range = extrema(p[2] for p in valid_points)
    z_range_grid = extrema(p[3] for p in valid_points)
    limits!(ax, x_range..., y_range..., z_range_grid...)

    # Pre-compute function values
    values = zeros(T, size(grid)...)
    @inbounds for i in eachindex(grid_points)
        point = grid_points[i]
        values[i] = any(isnan, point) ? NaN : f(point)
    end

    level_points = Observable(Point3f[])
    data_points = Observable(Point3f[])

    scatter!(ax, level_points, color = :blue, markersize = 2, label = "Level Set")

    if !isnothing(df)
        scatter!(
            ax,
            data_points,
            color = :darkorange,
            marker = :diamond,
            markersize = 20,
            label = "Data Points",
        )
    end

    function update_visualization(level::T) where {T<:AbstractFloat}
        level_data =
            prepare_level_set_data(grid, values, level, tolerance = params.point_tolerance)

        formatted_data = to_makie_format(level_data)

        new_points = Point3f[]
        if !isempty(formatted_data.xyz[1])
            for (x, y, z) in zip(formatted_data.xyz...)
                if !any(isnan, (x, y, z))
                    push!(new_points, Point3f(x, y, z))
                end
            end
        end
        level_points[] = new_points

        if !isnothing(df)
            visible_points = Point3f[]
            for row in eachrow(df)
                if !any(isnan, [row.x1, row.x2, row.x3, row.z]) &&
                   abs(row.z - level) ≤ params.point_tolerance
                    push!(visible_points, Point3f(row.x1, row.x2, row.x3))
                end
            end
            data_points[] = visible_points
        end
    end

    axislegend(ax, position = :rt)

    # Animation parameters
    total_frames = fps * duration
    θ = range(0, π, length = total_frames)
    levels = range(z_min, z_max, length = total_frames)

    record(fig, output_file, 1:total_frames; framerate = fps) do frame
        # Update camera position
        ax.azimuth[] = θ[frame]

        # Update level set
        update_visualization(levels[frame])
    end

    @info "Animation saved to: $output_file"
    return fig
end

# ═══════════════════════════════════════════════════════════════════════════════
# Pre-computed values overloads
#
# These methods accept an already-evaluated `values::Array{T}` instead of a
# callable `f`, avoiding redundant (and potentially expensive) sequential
# re-evaluation. Useful when values were computed via parallel grid evaluation
# (e.g., `evaluate_grid_threaded` from DynamicObjectives).
#
# Dispatch is unambiguous: Array{T} (these methods) vs callable (existing methods).
# ═══════════════════════════════════════════════════════════════════════════════

"""
    create_level_set_visualization(
        values::Array{T},
        grid::Array{SVector{3,T},3},
        df::Union{DataFrame,Nothing},
        z_range::Tuple{T,T},
        params::VisualizationParameters{T} = VisualizationParameters{T}()
    ) where {T<:AbstractFloat}

Create an interactive 3D level set visualization from pre-computed grid values.

Same as the `f`-based overload, but skips function evaluation — uses the provided
`values` array directly. The `values` array must have the same shape as `grid`.

`overlays` (keyword): point sets drawn while their value is within
`params.point_window` of the level, e.g. the critical points of several surrogates.
A vector of `(; df, label, marker = :circle, markersize = 16, color)`, where `df` has
columns x1, x2, x3, z on the scale of `values`, and optionally `index`, in which case
the points are coloured by Morse index. Sources and indices get a legend at
`fig[1, 2]`. `title` and `axis_labels` label the `Axis3`; `level_markersize` sets
the size of the level-set points (default 6 here, 2 in the animation).

`Inf` and `NaN` entries in `values` are treated as missing (excluded from level sets).
"""
function create_level_set_visualization(
    values::Array{T},
    grid::Array{SVector{3,T},3},
    df::Union{DataFrame,Nothing},
    z_range::Tuple{T,T},
    params::VisualizationParameters{T} = VisualizationParameters{T}();
    overlays = nothing,
    title::AbstractString = "",
    axis_labels = ("x₁", "x₂", "x₃"),
    level_markersize::Real = 6,
) where {T<:AbstractFloat}
    size(grid) == size(values) ||
        throw(DimensionMismatch("Grid shape $(size(grid)) != values shape $(size(values))"))

    grid_points = vec(grid)
    valid_points = filter(p -> !any(isnan, p), grid_points)
    isempty(valid_points) && throw(ArgumentError("Grid contains no valid points"))

    z_min, z_max = z_range
    (isnan(z_min) || isnan(z_max)) && throw(ArgumentError("Invalid z_range"))

    fig = Figure(size = params.fig_size)
    ax = Axis3(
        fig[1, 1],
        xlabel = axis_labels[1],
        ylabel = axis_labels[2],
        zlabel = axis_labels[3],
        title = title,
    )

    x_range = extrema(p[1] for p in valid_points)
    y_range = extrema(p[2] for p in valid_points)
    z_range_grid = extrema(p[3] for p in valid_points)

    limits!(ax, x_range..., y_range..., z_range_grid...)

    layers = overlays === nothing ? [] : _overlay_layers(ax, overlays)
    overlays === nothing || _overlay_legend!(fig[1, 2], overlays)

    level_slider =
        Slider(fig[2, 1], range = range(z_min, z_max, length = 1000), startvalue = z_min)

    Label(
        fig[3, 1],
        @lift(string("Level: ", round($(level_slider.value), digits = 3))),
        tellwidth = false,
    )

    level_points = Observable(Point3f[])
    data_points = Observable(Point3f[])

    scatter!(
        ax,
        level_points,
        color = :blue,
        markersize = level_markersize,
        alpha = 0.7,
        label = "Level Set",
    )

    if !isnothing(df)
        scatter!(
            ax,
            data_points,
            color = :darkorange,
            marker = :diamond,
            markersize = 30,
            label = "Data Points",
        )
    end

    function update_visualization(level::T) where {T<:AbstractFloat}
        try
            level_data = prepare_level_set_data(
                grid,
                values,
                level,
                tolerance = params.point_tolerance,
            )

            formatted_data = to_makie_format(level_data)

            new_points = Point3f[]
            if !isempty(formatted_data.xyz[1])
                for (x, y, z) in zip(formatted_data.xyz...)
                    if !any(isnan, (x, y, z))
                        push!(new_points, Point3f(x, y, z))
                    end
                end
            end
            level_points[] = new_points

            if !isnothing(df)
                visible_points = Point3f[]
                for row in eachrow(df)
                    if !any(isnan, [row.x1, row.x2, row.x3, row.z]) &&
                       abs(row.z - level) <= params.point_tolerance
                        push!(visible_points, Point3f(row.x1, row.x2, row.x3))
                    end
                end
                data_points[] = visible_points
            end
            _update_overlays!(layers, level, params.point_window)
        catch e
            @error "Error in visualization update" exception = e
            rethrow(e)
        end
    end

    on(level_slider.value) do level
        update_visualization(level)
    end

    update_visualization(z_min)
    overlays === nothing && axislegend(ax, position = :rt)

    return fig
end

"""
    create_level_set_animation(
        values::Array{T},
        grid::Array{SVector{3,T},3},
        df::Union{DataFrame,Nothing},
        z_range::Tuple{T,T},
        params::VisualizationParameters{T} = VisualizationParameters{T}();
        fps::Int = 30,
        duration::Int = 20,
        output_file::String = "level_set_animation.mp4"
    ) where {T<:AbstractFloat}

Create an animated level set visualization from pre-computed grid values.

Same as the `f`-based overload, but skips function evaluation — uses the provided
`values` array directly. The `values` array must have the same shape as `grid`.

`overlays` (keyword): point sets drawn while their value is within
`params.point_window` of the level, e.g. the critical points of several surrogates.
A vector of `(; df, label, marker = :circle, markersize = 16, color)`, where `df` has
columns x1, x2, x3, z on the scale of `values`, and optionally `index`, in which case
the points are coloured by Morse index. Sources and indices get a legend at
`fig[1, 2]`. `title` and `axis_labels` label the `Axis3`; `level_markersize` sets
the size of the level-set points (default 6 here, 2 in the animation).
"""
function create_level_set_animation(
    values::Array{T},
    grid::Array{SVector{3,T},3},
    df::Union{DataFrame,Nothing},
    z_range::Tuple{T,T},
    params::VisualizationParameters{T} = VisualizationParameters{T}();
    fps::Int = 30,
    duration::Int = 20,
    output_file::String = "level_set_animation.mp4",
    overlays = nothing,
    title::AbstractString = "Level Set Visualization",
    axis_labels = ("x₁", "x₂", "x₃"),
    level_markersize::Real = 2,
    level_label::Bool = overlays !== nothing,
    azimuth::Tuple{<:Real,<:Real} = (0, π),
    index_names = INDEX_NAMES,
    level_format = l -> string("level ", round(l, digits = 3)),
) where {T<:AbstractFloat}
    # `azimuth` is the camera's turn over the film (the default half turn passes edge-on at
    # π/2; a range strictly inside one quadrant, e.g. (1.05π, 1.45π), never does);
    # `index_names` relabels the Morse-index legend (e.g. valley floor / mountain pass /
    # hilltop); `level_format` writes the level caption, so a transformed value scale
    # (ranks, logs) can be captioned in the quantity's own units.
    size(grid) == size(values) ||
        throw(DimensionMismatch("Grid shape $(size(grid)) != values shape $(size(values))"))

    grid_points = vec(grid)
    valid_points = filter(p -> !any(isnan, p), grid_points)

    z_min, z_max = z_range

    fig = Figure(size = params.fig_size)
    ax = Axis3(
        fig[1, 1],
        title = title,
        xlabel = axis_labels[1],
        ylabel = axis_labels[2],
        zlabel = axis_labels[3],
    )

    # Set up initial ranges and limits
    x_range = extrema(p[1] for p in valid_points)
    y_range = extrema(p[2] for p in valid_points)
    z_range_grid = extrema(p[3] for p in valid_points)
    limits!(ax, x_range..., y_range..., z_range_grid...)

    layers = overlays === nothing ? [] : _overlay_layers(ax, overlays)
    overlays === nothing || _overlay_legend!(fig[1, 2], overlays; index_names = index_names)
    level_text = Observable("")
    level_label && Label(fig[2, 1], level_text, tellwidth = false)

    level_points = Observable(Point3f[])
    data_points = Observable(Point3f[])

    scatter!(
        ax,
        level_points,
        color = :blue,
        markersize = level_markersize,
        label = "Level Set",
    )

    if !isnothing(df)
        scatter!(
            ax,
            data_points,
            color = :darkorange,
            marker = :diamond,
            markersize = 20,
            label = "Data Points",
        )
    end

    function update_visualization(level::T) where {T<:AbstractFloat}
        level_data =
            prepare_level_set_data(grid, values, level, tolerance = params.point_tolerance)

        formatted_data = to_makie_format(level_data)

        new_points = Point3f[]
        if !isempty(formatted_data.xyz[1])
            for (x, y, z) in zip(formatted_data.xyz...)
                if !any(isnan, (x, y, z))
                    push!(new_points, Point3f(x, y, z))
                end
            end
        end
        level_points[] = new_points

        if !isnothing(df)
            visible_points = Point3f[]
            for row in eachrow(df)
                if !any(isnan, [row.x1, row.x2, row.x3, row.z]) &&
                   abs(row.z - level) <= params.point_tolerance
                    push!(visible_points, Point3f(row.x1, row.x2, row.x3))
                end
            end
            data_points[] = visible_points
        end
        _update_overlays!(layers, level, params.point_window)
        level_text[] = level_format(level)
    end

    overlays === nothing && axislegend(ax, position = :rt)

    # Animation parameters
    total_frames = fps * duration
    theta = range(azimuth[1], azimuth[2], length = total_frames)
    levels = range(z_min, z_max, length = total_frames)

    record(fig, output_file, 1:total_frames; framerate = fps) do frame
        # Update camera position
        ax.azimuth[] = theta[frame]

        # Update level set
        update_visualization(levels[frame])
    end

    @info "Animation saved to: $output_file"
    return fig
end

# ── sublevel-set animation (H₀ merge tree on the grid) ─────────────────────────

"""
    sublevel_barcode(values::AbstractArray{<:Real,3}) -> Vector{NamedTuple}

The H₀ barcode of the sublevel filtration of `values` on its grid (6-neighbourhood):
one `(; birth, death, argmin)` per connected component that ever appears, where
`birth` is its minimum value, `argmin` the linear index of that minimum, and `death`
the value at which it merges into a component with a lower minimum (elder rule);
the component of the global minimum has `death = Inf`. Non-finite values are skipped.
Sorted by birth.
"""
function sublevel_barcode(values::AbstractArray{<:Real,3})
    dims = size(values)
    N = length(values)
    LI = LinearIndices(dims)
    CI = CartesianIndices(dims)
    parent = zeros(Int, N)
    offs = (
        CartesianIndex(1, 0, 0),
        CartesianIndex(-1, 0, 0),
        CartesianIndex(0, 1, 0),
        CartesianIndex(0, -1, 0),
        CartesianIndex(0, 0, 1),
        CartesianIndex(0, 0, -1),
    )
    function findroot(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    death = Dict{Int,Float64}()
    for v in sortperm(vec(values))
        isfinite(values[v]) || continue
        parent[v] = v
        for o in offs
            u = CI[v] + o
            checkbounds(Bool, values, u) || continue
            ui = LI[u]
            parent[ui] == 0 && continue
            ru, rv = findroot(ui), findroot(v)
            ru == rv && continue
            # roots are the components' minima: the higher one dies here
            old, young = values[ru] <= values[rv] ? (ru, rv) : (rv, ru)
            parent[young] = old
            young == v || (death[young] = Float64(values[v]))
        end
    end
    roots = [v for v in 1:N if parent[v] == v]      # surviving components (normally one)
    bars = [(; birth = Float64(values[k]), death = d, argmin = k) for (k, d) in death]
    append!(bars, [(; birth = Float64(values[k]), death = Inf, argmin = k) for k in roots])
    return sort(bars; by = b -> b.birth)
end

const COMPONENT_PALETTE = (
    :royalblue,
    :darkorange,
    :forestgreen,
    :crimson,
    :purple,
    :saddlebrown,
    :deeppink,
    :gray40,
    :olive,
    :teal,
    :gold,
    :navy,
)

"""
    create_sublevel_animation(values, grid, z_range, params = VisualizationParameters();
        overlays = nothing, fps = 30, duration = 20, output_file = "sublevel_animation.mp4",
        title = "Sublevel sets", axis_labels = ("x₁", "x₂", "x₃"), level_markersize = 2,
        rotate = true) -> Figure

Animate the sublevel sets {values ≤ c} as c sweeps `z_range`. The sets are grown by
union–find over the grid's 6-neighbourhood (the H₀ merge tree), so each connected
component is known at every level. The boundary shell
(c − `params.point_tolerance` ≤ values ≤ c) is drawn in the colour of its
component. A component's colour is fixed by the rank of its minimum, and when two
components merge the one with the lower minimum keeps its colour (elder rule), so a
merge shows as one colour swallowing another. The number of components is shown
under the axis and, with `barcode = true`, the H₀ barcode (`sublevel_barcode`, bars
shorter than `min_persistence` dropped) with a cursor at c.
`overlays` are as in `create_level_set_visualization` (drawn within
`params.point_window` of c): saddles of index 1 light up at the merge levels.
"""
function create_sublevel_animation(
    values::Array{T,3},
    grid::Array{SVector{3,T},3},
    z_range::Tuple{T,T},
    params::VisualizationParameters{T} = VisualizationParameters{T}();
    overlays = nothing,
    fps::Int = 30,
    duration::Int = 20,
    output_file::String = "sublevel_animation.mp4",
    title::AbstractString = "Sublevel sets",
    axis_labels = ("x₁", "x₂", "x₃"),
    level_markersize::Real = 2,
    rotate::Bool = true,
    barcode::Bool = true,
    min_persistence::Real = 0,
) where {T<:AbstractFloat}
    size(grid) == size(values) ||
        throw(DimensionMismatch("Grid shape $(size(grid)) != values shape $(size(values))"))
    z_min, z_max = z_range
    z_min < z_max || throw(ArgumentError("z_range must be increasing, got $z_range"))
    dims = size(values)
    N = length(values)
    LI = LinearIndices(dims)
    CI = CartesianIndices(dims)
    order = sortperm(vec(values))
    parent = zeros(Int, N)                  # 0 = not yet in the sublevel set
    birth = zeros(Int, N)                   # birth rank of the component whose root this is
    nborn = Ref(0)
    ncomp = Ref(0)
    function findroot(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    offs = (
        CartesianIndex(1, 0, 0),
        CartesianIndex(-1, 0, 0),
        CartesianIndex(0, 1, 0),
        CartesianIndex(0, -1, 0),
        CartesianIndex(0, 0, 1),
        CartesianIndex(0, 0, -1),
    )
    ptr = Ref(1)
    function grow_to!(c)
        while ptr[] <= N && values[order[ptr[]]] <= c
            v = order[ptr[]]
            isfinite(values[v]) || (ptr[] += 1; continue)
            parent[v] = v
            nborn[] += 1
            birth[v] = nborn[]
            ncomp[] += 1
            for o in offs
                u = CI[v] + o
                checkbounds(Bool, values, u) || continue
                ui = LI[u]
                parent[ui] == 0 && continue
                ru, rv = findroot(ui), findroot(v)
                ru == rv && continue
                # elder rule: the component born first (lower minimum) survives
                old, young = birth[ru] < birth[rv] ? (ru, rv) : (rv, ru)
                parent[young] = old
                ncomp[] -= 1
            end
            ptr[] += 1
        end
    end

    # a component is named by its minimum (the union-find root, elder rule); its colour
    # is its rank in the barcode, so shell and bar share a colour
    allbars = sublevel_barcode(values)
    rank = Dict(b.argmin => k for (k, b) in enumerate(allbars))
    grid_points = vec(grid)
    fig = Figure(size = params.fig_size)
    ax = Axis3(
        fig[1, 1],
        title = title,
        xlabel = axis_labels[1],
        ylabel = axis_labels[2],
        zlabel = axis_labels[3],
    )
    limits!(
        ax,
        extrema(p[1] for p in grid_points)...,
        extrema(p[2] for p in grid_points)...,
        extrema(p[3] for p in grid_points)...,
    )
    P = length(COMPONENT_PALETTE)
    slots = [Observable(Point3f[]) for _ in 1:P]
    for (k, obs) in enumerate(slots)
        scatter!(ax, obs; color = COMPONENT_PALETTE[k], markersize = level_markersize)
    end
    layers = overlays === nothing ? [] : _overlay_layers(ax, overlays)
    overlays === nothing || _overlay_legend!(fig[1, 2], overlays)
    info = Observable("")
    Label(fig[2, 1], info, tellwidth = false)

    # barcode: one bar per component of {f ≤ c} over the animated range, birth → merge
    # (elder rule), bars of persistence < min_persistence dropped; a cursor at c
    cursor = Observable(z_min)
    if barcode
        bars = [
            b for b in allbars if
            b.birth <= z_max && min(b.death, z_max) - b.birth >= min_persistence
        ]
        bx = Axis(
            fig[3, 1],
            xlabel = "level",
            ylabel = "component",
            height = 140,
            yticksvisible = false,
            yticklabelsvisible = false,
        )
        xlims!(bx, z_min, z_max)
        for (k, b) in enumerate(bars)
            lines!(
                bx,
                [b.birth, min(b.death, z_max)],
                [k, k];
                color = COMPONENT_PALETTE[mod1(rank[b.argmin], length(COMPONENT_PALETTE))],
                linewidth = 4,
            )
        end
        isempty(bars) || ylims!(bx, 0, length(bars) + 1)
        vlines!(bx, cursor; color = :black, linewidth = 1)
    end

    tol = params.point_tolerance
    function draw!(c)
        grow_to!(c)
        buckets = [Point3f[] for _ in 1:P]
        for v in eachindex(values)
            parent[v] == 0 && continue
            values[v] >= c - tol || continue
            r = findroot(v)
            p = grid_points[v]
            push!(buckets[mod1(rank[r], P)], Point3f(p[1], p[2], p[3]))
        end
        for k in 1:P
            slots[k][] = buckets[k]
        end
        _update_overlays!(layers, c, params.point_window)
        info[] = string(
            "level ",
            round(c, digits = 3),
            "    components of {f ≤ level}: ",
            ncomp[],
        )
        cursor[] = c
    end

    total_frames = fps * duration
    theta = range(0, pi, length = total_frames)
    levels = range(z_min, z_max, length = total_frames)
    record(fig, output_file, 1:total_frames; framerate = fps) do frame
        rotate && (ax.azimuth[] = theta[frame])
        draw!(levels[frame])
    end
    @info "Sublevel animation saved to: $output_file"
    return fig
end

# ── parameter-sweep animation (one grid per parameter value, fixed level) ─────

"""
    create_level_set_sweep_animation(frames, grid, level, params = VisualizationParameters();
        frame_overlays = nothing, frame_labels = nothing, substeps = 8, fps = 30,
        output_file = "level_set_sweep.mp4", title = "", axis_labels = ("x₁", "x₂", "x₃"),
        level_markersize = 2, rotate = true) -> Figure

Animate a family of functions on one grid, e.g. slices f(·; t) of a higher-dimensional
objective at parameter values t₁ < … < t_K, at the fixed `level`. `frames[j]` holds the
values of the j-th member on `grid`. Between consecutive members the values are
interpolated linearly in `substeps` steps, so the level shell moves continuously.
`frame_overlays[j]` (overlay vectors as in `create_level_set_visualization`, the same
sources in the same order in every frame) is drawn
while frame j is the nearest member, so critical points appear and vanish where the
family has a fold. `params.point_window` decides which overlay points are shown
(use a large window to show them all). `frame_labels[j]` is shown under the axis.
"""
function create_level_set_sweep_animation(
    frames::AbstractVector{<:Array{T,3}},
    grid::Array{SVector{3,T},3},
    level::Real,
    params::VisualizationParameters{T} = VisualizationParameters{T}();
    frame_overlays = nothing,
    frame_labels = nothing,
    substeps::Int = 8,
    fps::Int = 30,
    output_file::String = "level_set_sweep.mp4",
    title::AbstractString = "",
    axis_labels = ("x₁", "x₂", "x₃"),
    level_markersize::Real = 2,
    rotate::Bool = true,
) where {T<:AbstractFloat}
    K = length(frames)
    K >= 1 || throw(ArgumentError("frames is empty"))
    all(size(F) == size(grid) for F in frames) ||
        throw(DimensionMismatch("every frame must have the grid's shape $(size(grid))"))
    frame_overlays === nothing ||
        length(frame_overlays) == K ||
        throw(ArgumentError("frame_overlays needs one entry per frame ($K)"))
    frame_labels === nothing ||
        length(frame_labels) == K ||
        throw(ArgumentError("frame_labels needs one entry per frame ($K)"))
    substeps >= 1 || throw(ArgumentError("substeps must be ≥ 1"))

    grid_points = vec(grid)
    fig = Figure(size = params.fig_size)
    ax = Axis3(
        fig[1, 1],
        title = title,
        xlabel = axis_labels[1],
        ylabel = axis_labels[2],
        zlabel = axis_labels[3],
    )
    limits!(
        ax,
        extrema(p[1] for p in grid_points)...,
        extrema(p[2] for p in grid_points)...,
        extrema(p[3] for p in grid_points)...,
    )
    shell = Observable(Point3f[])
    scatter!(ax, shell; color = :blue, markersize = level_markersize)
    # one set of overlay layers per frame, all built up front; the inactive ones are emptied
    layer_sets =
        frame_overlays === nothing ? [] : [_overlay_layers(ax, ov) for ov in frame_overlays]
    if frame_overlays !== nothing
        # legend: the sources of the first frame (every frame has the same sources),
        # with the Morse indices pooled over all frames
        first_ov = frame_overlays[1]
        pooled = [
            merge(
                first_ov[i],
                (;
                    df = reduce(
                        vcat,
                        [ov[i].df for ov in frame_overlays if length(ov) >= i],
                    )
                ),
            ) for i in eachindex(first_ov)
        ]
        _overlay_legend!(fig[1, 2], pooled)
    end
    info = Observable("")
    Label(fig[2, 1], info, tellwidth = false)

    tol = params.point_tolerance
    steps = [
        (j, a) for j in 1:max(K - 1, 1) for a in range(0, 1; length = substeps + 1)[1:end-1]
    ]
    K > 1 && push!(steps, (K - 1, 1.0))
    theta = range(0, pi, length = length(steps))
    record(fig, output_file, eachindex(steps); framerate = fps) do s
        j, a = steps[s]
        A = frames[j]
        B = K > 1 ? frames[j+1] : frames[j]
        pts = Point3f[]
        for v in eachindex(A)
            x = (1 - a) * A[v] + a * B[v]
            (isfinite(x) && abs(x - level) <= tol) || continue
            p = grid_points[v]
            push!(pts, Point3f(p[1], p[2], p[3]))
        end
        shell[] = pts
        near = K > 1 && a >= 0.5 ? j + 1 : j
        for (m, layers) in enumerate(layer_sets)
            _update_overlays!(layers, level, m == near ? params.point_window : -one(T))
        end
        frame_labels === nothing || (info[] = frame_labels[near])
        rotate && (ax.azimuth[] = theta[s])
    end
    @info "Sweep animation saved to: $output_file"
    return fig
end
