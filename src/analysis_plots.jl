using LinearAlgebra
using DataFrames
using Statistics

"""
Analyze convergence distances for a DataFrame of critical points.
Returns statistics about the distances between points.
"""
function analyze_convergence_distances(df::DataFrame)
    # Get dimension from column names
    dim = count(col -> startswith(string(col), "x"), names(df))

    # Calculate minimum distance for each point to any other point
    n_points = nrow(df)
    min_distances = Float64[]

    for i in 1:n_points
        point_i = [df[i, Symbol("x$j")] for j in 1:dim]
        min_dist = Inf

        for j in 1:n_points
            if i != j
                point_j = [df[j, Symbol("x$k")] for k in 1:dim]
                dist = norm(point_i - point_j)
                min_dist = min(min_dist, dist)
            end
        end

        if min_dist < Inf
            push!(min_distances, min_dist)
        end
    end

    # Return statistics
    if isempty(min_distances)
        return (maximum = 0.0, average = 0.0, minimum = 0.0)
    else
        return (
            maximum = maximum(min_distances),
            average = mean(min_distances),
            minimum = minimum(min_distances),
        )
    end
end

"""
Analyze distances between captured points and reference points.
Returns statistics about minimum distances from each point in df to closest point in df_check.
"""
function analyze_captured_distances(df::DataFrame, df_check::DataFrame)
    # Get dimension from column names
    dim = count(col -> startswith(string(col), "x"), names(df))

    # Calculate minimum distance for each point in df to any point in df_check
    min_distances = Float64[]

    for i in 1:nrow(df)
        point = [df[i, Symbol("x$k")] for k in 1:dim]
        min_dist = Inf

        for j in 1:nrow(df_check)
            check_point = [df_check[j, Symbol("x$k")] for k in 1:dim]
            dist = norm(point - check_point)
            min_dist = min(min_dist, dist)
        end

        push!(min_distances, min_dist)
    end

    # Return statistics
    if isempty(min_distances)
        return (maximum = 0.0, average = 0.0, minimum = 0.0)
    else
        return (
            maximum = maximum(min_distances),
            average = mean(min_distances),
            minimum = minimum(min_distances),
        )
    end
end

"""
Plot the discrete L2-norm approximation error attained by the polynomial approximant.
"""
function plot_discrete_l2(
    results,
    start_degree::Int,
    end_degree::Int,
    step::Int;
    fig_size::Tuple{Int,Int} = (700, 500),
)
    # Filter to only include degrees that succeeded
    all_degrees = start_degree:step:end_degree
    degrees = filter(d -> haskey(results, d), all_degrees)

    if isempty(degrees)
        error("No successful results found for degrees $start_degree:$step:$end_degree")
    end

    l2_norms = Float64[]

    # Extract L2 norms for each degree
    for d in degrees
        push!(l2_norms, results[d].discrete_l2)
    end

    # Create figure
    fig = Figure(size = fig_size, fontsize = 14)

    ax = Axis(
        fig[1, 1],
        xlabel = "Polynomial Degree",
        ylabel = "L² Approximation Error",
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xticks = degrees,
    )

    # Plot the curve with points at each degree
    scatterlines!(
        ax,
        degrees,
        l2_norms,
        color = :darkblue,
        markersize = 10,
        linewidth = 2.5,
        label = "L² Norm",
    )

    axislegend(ax, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

    return fig
end

"""
Plot summary of convergence distances for a range of degrees --> for each captured "x", compute the distance to "y", the optimized point.
"""
function plot_convergence_analysis(
    results,
    start_degree::Int,
    end_degree::Int,
    step::Int;
    show_legend::Bool = true,
    fig_size::Tuple{Int,Int} = (700, 500),
)
    # Filter to only include degrees that succeeded
    all_degrees = start_degree:step:end_degree
    degrees = filter(d -> haskey(results, d), all_degrees)

    if isempty(degrees)
        error("No successful results found for degrees $start_degree:$step:$end_degree")
    end

    max_distances = Float64[]
    avg_distances = Float64[]

    for d in degrees
        df = results[d].df
        stats = analyze_convergence_distances(df)
        push!(max_distances, stats.maximum)
        push!(avg_distances, stats.average)
    end

    fig = Figure(size = fig_size, fontsize = 14)

    ax = Axis(
        fig[1, 1],
        xlabel = "Polynomial Degree",
        ylabel = "Distance to Nearest Critical Point",
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xticks = degrees,
    )

    scatterlines!(
        ax,
        degrees,
        max_distances,
        label = "Maximum",
        color = :crimson,
        markersize = 10,
        linewidth = 2.5,
    )
    scatterlines!(
        ax,
        degrees,
        avg_distances,
        label = "Average",
        color = :steelblue,
        markersize = 10,
        linewidth = 2.5,
    )

    axislegend(ax, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

    return fig
end

"""
    cairo_plot_polyapprox_levelset(pol, TR, df, df_min; kwargs...)

Plot polynomial approximation level sets with critical points overlaid.

Handles per-coordinate scaling factors. Creates a CairoMakie contour plot
with scattered critical points (near/far, captured/uncaptured).

# Keyword arguments
- `figure_size::Tuple{Int,Int} = (1000, 600)`: figure dimensions
- `z_limits::Union{Nothing, Tuple{Float64,Float64}} = nothing`: contour z range (auto if nothing)
- `chebyshev_levels::Bool = false`: use Chebyshev-spaced contour levels
- `num_levels::Int = 30`: number of contour levels
- `show_captured::Bool = true`: whether to show captured minima
"""
function cairo_plot_polyapprox_levelset(
    pol::AbstractPolynomialData,
    TR::AbstractProblemInput,
    df::DataFrame,
    df_min::DataFrame;
    figure_size::Tuple{Int,Int} = (1000, 600),
    z_limits::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    chebyshev_levels::Bool = false,
    num_levels::Int = 30,
    show_captured::Bool = true,
    colormap::Symbol = :viridis,
    title::AbstractString = "",
    title_fontsize::Real = 14,
    xlabel::AbstractString = "x₁",
    ylabel::AbstractString = "x₂",
    colorbar_label::AbstractString = "f(x)",
    show_colorbar::Bool = true,
    show_legend::Bool = true,
    legend_below::Bool = false,
    label_far::AbstractString = "Far",
    label_near::AbstractString = "Near",
    label_captured::AbstractString = "Captured",
    label_uncaptured::AbstractString = "Uncaptured",
)
    # Type-stable coordinate transformation using multiple dispatch
    coords = transform_coordinates(pol.scale_factor, pol.grid, TR.center)

    z_coords = pol.z

    if size(coords)[2] == 2
        fig = CairoMakie.Figure(size = figure_size)
        ax = CairoMakie.Axis(
            fig[1, 1];
            title = title,
            titlesize = title_fontsize,
            xlabel = xlabel,
            ylabel = ylabel,
        )

        # Calculate z_limits if not provided (filter non-finite values from ODE solver failures)
        if isnothing(z_limits)
            z_values = Float64[]
            append!(z_values, df.z)
            append!(z_values, df_min.value)
            finite_z = filter(isfinite, z_values)
            if isempty(finite_z)
                error(
                    "No finite z-values among critical points — cannot determine contour z_limits",
                )
            end
            z_limits = (minimum(finite_z), maximum(finite_z))
        end

        # Calculate levels (must be finite and sorted for Makie contourf)
        levels = if chebyshev_levels
            k = collect(0:(num_levels-1))
            cheb_nodes = -cos.((2k .+ 1) .* π ./ (2 * num_levels))
            z_min, z_max = z_limits
            lvls = (z_max - z_min) ./ 2 .* cheb_nodes .+ (z_max + z_min) ./ 2
            sort!(filter!(isfinite, lvls))
        else
            num_levels
        end

        # Prepare contour data
        x_unique = sort(unique(coords[:, 1]))
        y_unique = sort(unique(coords[:, 2]))
        Z = fill(NaN, (length(y_unique), length(x_unique)))

        for (idx, (x, y, z)) in enumerate(zip(coords[:, 1], coords[:, 2], z_coords))
            i = findlast(≈(y), y_unique)
            j = findlast(≈(x), x_unique)
            if !isnothing(i) && !isnothing(j)
                Z[j, i] = isfinite(z) ? z : NaN  # non-finite → NaN (empty cell in contourf)
            end
        end

        # Create contour plot
        chosen_colormap = colormap
        cf = CairoMakie.contourf!(
            ax,
            x_unique,
            y_unique,
            Z,
            colormap = chosen_colormap,
            levels = levels,
        )

        # Initialize empty array for legend entries
        legend_entries = []

        # Plot and add legend entries for all point types
        if :close in propertynames(df)
            # Far points
            not_close_idx = .!df.close
            if any(not_close_idx)
                CairoMakie.scatter!(
                    ax,
                    df.x1[not_close_idx],
                    df.x2[not_close_idx],
                    markersize = 10,
                    color = :white,
                    strokecolor = :black,
                    strokewidth = 1,
                    label = label_far,
                )
                push!(legend_entries, label_far)
            end

            # Near points
            close_idx = df.close
            if any(close_idx)
                CairoMakie.scatter!(
                    ax,
                    df.x1[close_idx],
                    df.x2[close_idx],
                    markersize = 10,
                    color = :green,
                    strokecolor = :black,
                    strokewidth = 1,
                    label = label_near,
                )
                push!(legend_entries, label_near)
            end
        else
            # All points if no close/far distinction
            CairoMakie.scatter!(
                ax,
                df.x1,
                df.x2,
                markersize = 2,
                color = :orange,
                label = "All points",
            )
            push!(legend_entries, "All points")
        end

        # Uncaptured points
        if !isempty(df_min)
            uncaptured_idx = .!df_min.captured
            captured_idx = df_min.captured

            if any(uncaptured_idx)
                CairoMakie.scatter!(
                    ax,
                    df_min.x1[uncaptured_idx],
                    df_min.x2[uncaptured_idx],
                    markersize = 15,
                    marker = :diamond,
                    color = :red,
                    label = label_uncaptured,
                )
                push!(legend_entries, label_uncaptured)
            end

            # Only show captured points if show_captured is true
            if show_captured && any(captured_idx)
                CairoMakie.scatter!(
                    ax,
                    df_min.x1[captured_idx],
                    df_min.x2[captured_idx],
                    markersize = 15,
                    marker = :diamond,
                    color = :blue,
                    label = label_captured,
                )
                push!(legend_entries, label_captured)
            end
        end

        if show_colorbar
            CairoMakie.Colorbar(fig[1, 2], cf; label = colorbar_label, width = 14)
        end
        if show_legend && !isempty(legend_entries)
            if legend_below
                CairoMakie.Legend(
                    fig[2, 1:(show_colorbar ? 2 : 1)],
                    ax;
                    orientation = :horizontal,
                    framevisible = false,
                    labelsize = 12,
                    patchsize = (18, 18),
                )
            else
                CairoMakie.axislegend(
                    ax;
                    position = :rt,
                    framevisible = true,
                    backgroundcolor = (:white, 0.85),
                    labelsize = 11,
                )
            end
        end
        return fig
    end
end

"""
Plot the outputs of`analyze_converged_points` function.
"""
function plot_distance_statistics(
    stats::Dict{String,Any};
    show_legend::Bool = true,
    fig_size::Tuple{Int,Int} = (700, 500),
)
    fig = Figure(size = fig_size, fontsize = 14)

    degrees = stats["degrees"]

    ax = Axis(
        fig[1, 1],
        xlabel = "Polynomial Degree",
        ylabel = "Distance",
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xticks = degrees,
    )

    # Plot maximum and average distances
    scatterlines!(
        ax,
        degrees,
        stats["max_distances"],
        label = "Maximum",
        color = :crimson,
        markersize = 10,
        linewidth = 2.5,
    )
    scatterlines!(
        ax,
        degrees,
        stats["avg_distances"],
        label = "Average",
        color = :steelblue,
        markersize = 10,
        linewidth = 2.5,
    )

    axislegend(ax, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

    return fig
end

"""
    plot_convergence_captured(results, df_check, start_degree, end_degree, step; kwargs...) -> Figure

Distance from the critical points found at each degree to the nearest reference point
of `df_check` (e.g. known minimizers), as maximum and average against the degree.
`results` is a `Dict` keyed by degree whose values carry the critical-point frame in
`.df`; degrees missing from it are skipped.
"""
function plot_convergence_captured(
    results,
    df_check,
    start_degree::Int,
    end_degree::Int,
    step::Int;
    show_legend::Bool = true,
    fig_size::Tuple{Int,Int} = (600, 400),
)
    # Filter to only include degrees that succeeded
    all_degrees = start_degree:step:end_degree
    degrees = filter(d -> haskey(results, d), all_degrees)

    if isempty(degrees)
        error("No successful results found for degrees $start_degree:$step:$end_degree")
    end

    max_distances = Float64[]
    avg_distances = Float64[]

    for d in degrees
        x_cols = [col for col in names(results[d].df_min) if startswith(string(col), "x")]
        df = results[d].df_min[:, x_cols]

        stats = analyze_captured_distances(df, df_check)
        push!(max_distances, stats.maximum)
        push!(avg_distances, stats.average)
    end

    fig = Figure(size = fig_size)

    ax = Axis(
        fig[1, 1],
        # title="",
        xlabel = "Degree",
        ylabel = "",
    )

    scatterlines!(ax, degrees, max_distances, label = "Maximum", color = :red)
    scatterlines!(ax, degrees, avg_distances, label = "Average", color = :blue)

    # Legend removed per user request

    return fig
end

"""
Plot distance from theoretical minimizers to nearest computed critical points.

For each theoretical minimizer, computes the distance to the nearest critical point
found by the polynomial method. Shows both average and maximum distances by degree.

# Arguments
- `results`: Dictionary mapping degrees to analysis results
- `df_theoretical_mins`: DataFrame containing theoretical minimizer coordinates (x1, x2, ...)
- `start_degree`: Minimum polynomial degree
- `end_degree`: Maximum polynomial degree
- `step`: Step size between degrees (default: 1)
- `show_legend`: Whether to display legend (default: true)

# Returns
- `fig`: Makie Figure object
"""
function plot_theoretical_minimizer_distances(
    results,
    df_theoretical_mins::DataFrame,
    start_degree::Int,
    end_degree::Int,
    step::Int = 1;
    show_legend::Bool = true,
    fig_size::Tuple{Int,Int} = (700, 500),
)
    # Filter to only include degrees that succeeded
    all_degrees = start_degree:step:end_degree
    degrees = filter(d -> haskey(results, d), all_degrees)

    if isempty(degrees)
        error("No successful results found for degrees $start_degree:$step:$end_degree")
    end

    # Determine dimension from theoretical minimizers
    dim = count(col -> startswith(string(col), "x"), names(df_theoretical_mins))

    avg_distances = Float64[]
    max_distances = Float64[]

    for d in degrees
        df_computed = results[d].df  # All computed critical points
        distances = Float64[]

        # For each theoretical minimizer, find distance to nearest computed point
        for i in 1:nrow(df_theoretical_mins)
            theo_pt = [df_theoretical_mins[i, Symbol("x$j")] for j in 1:dim]
            min_dist = Inf

            # Find nearest computed critical point
            for j in 1:nrow(df_computed)
                computed_pt = [df_computed[j, Symbol("x$j")] for j in 1:dim]
                dist = norm(theo_pt - computed_pt)
                min_dist = min(min_dist, dist)
            end

            push!(distances, min_dist)
        end

        # Compute statistics for this degree
        push!(avg_distances, mean(distances))
        push!(max_distances, maximum(distances))
    end

    # Create figure
    fig = Figure(size = fig_size, fontsize = 14)

    ax = Axis(
        fig[1, 1],
        xlabel = "Polynomial Degree",
        ylabel = "Distance to Theoretical Minimizers",
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
        xticks = degrees,
    )

    scatterlines!(
        ax,
        degrees,
        avg_distances,
        label = "Average",
        color = :steelblue,
        markersize = 10,
        linewidth = 2.5,
    )
    scatterlines!(
        ax,
        degrees,
        max_distances,
        label = "Maximum",
        color = :crimson,
        markersize = 10,
        linewidth = 2.5,
    )

    axislegend(ax, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

    return fig
end

"""
    plot_l2_convergence_trials(df::DataFrame; kwargs...)

Plot L2 approximation error vs polynomial degree with individual trial paths.

# Arguments
- `df::DataFrame`: Must have columns `:degree`, `:l2_error`, `:seed`, `:method`

# Keyword Arguments
- `domain_filter::Union{Nothing, Float64}=nothing`: Filter to specific domain
- `title::String="L2 Approximation Error by Degree"`: Plot title
- `fig_size::Tuple{Int,Int}=(800, 500)`: Figure size
- `show_median::Bool=true`: Overlay median line per method
- `trial_alpha::Float64=0.3`: Alpha for individual trial lines
- `colors::Tuple=(Makie.wong_colors()[1], Makie.wong_colors()[2])`: (standard, log) colors

# Returns
- `Figure`: Makie figure object
"""
function plot_l2_convergence_trials(
    df::DataFrame;
    domain_filter::Union{Nothing,Float64} = nothing,
    title::String = "L2 Approximation Error by Degree",
    fig_size::Tuple{Int,Int} = (800, 500),
    show_median::Bool = true,
    trial_alpha::Float64 = 0.3,
    colors = (Makie.wong_colors()[1], Makie.wong_colors()[2]),
)
    # Filter by domain if specified
    data = isnothing(domain_filter) ? df : filter(r -> r.domain == domain_filter, df)

    fig = Figure(size = fig_size, fontsize = 14)
    ax = Axis(
        fig[1, 1],
        title = title,
        xlabel = "Polynomial Degree",
        ylabel = "L2 Error",
        yscale = log10,
        xgridvisible = true,
        ygridvisible = true,
        xgridstyle = :dash,
        ygridstyle = :dash,
    )

    method_colors = Dict("standard" => colors[1], "log" => colors[2])

    # Plot individual trial paths
    for method in ["standard", "log"]
        method_data = filter(r -> r.method == method, data)
        isempty(method_data) && continue

        color = method_colors[method]
        seeds = unique(method_data.seed)

        for seed in seeds
            trial = filter(r -> r.seed == seed, method_data)
            sort!(trial, :degree)
            lines!(
                ax,
                trial.degree,
                trial.l2_error,
                color = (color, trial_alpha),
                linewidth = 1.5,
            )
        end

        # Median line per method
        if show_median
            degrees = sort(unique(method_data.degree))
            medians =
                [median(filter(r -> r.degree == d, method_data).l2_error) for d in degrees]
            scatterlines!(
                ax,
                degrees,
                medians,
                color = color,
                linewidth = 3,
                markersize = 10,
                label = method,
            )
        end
    end

    axislegend(ax, position = :rt, framevisible = true, bgcolor = (:white, 0.9))

    return fig
end
