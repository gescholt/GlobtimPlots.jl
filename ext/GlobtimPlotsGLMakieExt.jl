module GlobtimPlotsGLMakieExt

# GLMakie extension for GlobtimPlots.
#
# All plotting functions live in GlobtimPlots/src using the abstract Makie API and
# work with whichever backend is active. This extension adds the GLMakie-specific
# `plot_polyapprox_3d` method: the backend-agnostic surface plus an optional
# recorded rotation animation (`rotate = true`, needs GLMakie.record).

using GlobtimPlots
using GLMakie
using DataFrames

import GlobtimPlots: AbstractPolynomialData, AbstractProblemInput

"""
    plot_polyapprox_3d(pol, TR, df, df_min; kwargs...)

3D polynomial approximation surface with critical points (GLMakie backend).
With `rotate = true` a 240-frame rotation is recorded to `filename` via `GLMakie.record`.

# Returns
- `fig`: The GLMakie figure object
"""
function GlobtimPlots.plot_polyapprox_3d(
    pol::AbstractPolynomialData,
    TR::AbstractProblemInput,
    df::DataFrame,
    df_min::DataFrame;
    figure_size::Tuple{Int,Int} = (1000, 800),
    z_limits::Union{Nothing,Tuple{Float64,Float64}} = nothing,
    show_captured::Bool = true,
    alpha_surface::Float64 = 0.7,
    rotate::Bool = false,
    filename::String = "function_3d_rotation.mp4",
    fade::Bool = false,
    z_cut = 0.25,
    color_by::Symbol = :proximity,
)
    # Delegate to the backend-agnostic implementation in GlobtimPlots
    fig = GlobtimPlots._plot_polyapprox_3d_impl(
        pol,
        TR,
        df,
        df_min;
        figure_size = figure_size,
        z_limits = z_limits,
        show_captured = show_captured,
        alpha_surface = alpha_surface,
        fade = fade,
        z_cut = z_cut,
        color_by = color_by,
    )

    if isnothing(fig)
        return nothing
    end

    # GLMakie-specific: optional rotation animation (requires GLMakie.record)
    if rotate
        ax = content(fig[1, 1])
        @info "Recording animation to $(filename)..."
        GLMakie.record(fig, filename, 1:240; framerate = 30) do frame
            ax.azimuth[] = 3π / 4 + 2π * frame / 240
            ax.elevation[] = π / 6 + π / 12 * sin(2π * frame / 240)
        end
        @info "Animation saved to $(filename)!"
    end

    display(fig)
    return fig
end

end
