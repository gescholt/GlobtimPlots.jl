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
function _overlay_legend!(pos, overlays)
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
        push!(labels, INDEX_NAMES[k+1])
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
) where {T<:AbstractFloat}
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
    overlays === nothing || _overlay_legend!(fig[1, 2], overlays)
    level_text = Observable("")
    level_label && Label(fig[2, 1], level_text, tellwidth = false)

    level_points = Observable(Point3f[])
    data_points = Observable(Point3f[])

    scatter!(ax, level_points, color = :blue, markersize = level_markersize, label = "Level Set")

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
        level_text[] = string("level ", round(level, digits = 3))
    end

    overlays === nothing && axislegend(ax, position = :rt)

    # Animation parameters
    total_frames = fps * duration
    theta = range(0, pi, length = total_frames)
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
