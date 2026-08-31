# =============================================================================
# paper_theme.jl — shared figure style for the CertifiedParameterEstimation
# manuscript (plain `article` class, no geometry package: \linewidth = 345 pt).
#
# The historical failure mode this fixes: figures authored at 1060 Makie units
# and included at width=\linewidth get scaled by 345/795 ≈ 0.43, so fontsize 18
# lands at ~5.9 pt on the page — and every script picked its own size, so
# effective font sizes varied ~1.5x across the paper. The convention here is to
# author at TRUE PAGE SIZE: CairoMakie saves PDFs at 0.75 pt per Makie unit, so
# a figure `paper_figsize(1.0)` = 460 units wide is exactly 345 pt and is
# included 1:1 with `\includegraphics[width=\linewidth]` — fonts on the page
# are exactly the theme's sizes (9 pt body, 8 pt ticks).
#
# Usage in a figure script:
#     using GlobtimPlots   # or: include this file
#     with_theme(paper_theme()) do
#         fig = Figure(size = paper_figsize(1.0; aspect = 0.55))
#         ...
#         save("fig.pdf", fig)                 # vector, true size
#         save("fig.png", fig; px_per_unit = 4)  # preview/raster twin
#     end
# =============================================================================

export paper_theme, paper_figsize, PAPER_LINEWIDTH_UNITS, PAPER_PALETTE

"\\linewidth of the manuscript (plain 10pt article: 345 pt) in Makie units (0.75 pt/unit)."
const PAPER_LINEWIDTH_UNITS = 460.0

"House palette: teal / red / orange / grey — matches the existing dcert/valley figures."
const PAPER_PALETTE = ["#1b7f8e", "#c8503c", "#e08214", "#5e6a71", "#7b3294", "#008837"]

"""
    paper_figsize(fraction = 1.0; aspect = 0.55) -> (w, h)

Figure size in Makie units for a figure included at `fraction * \\linewidth`.
`aspect` is height/width. Author at this size and include at the SAME fraction
so page font sizes equal the theme's sizes.
"""
function paper_figsize(fraction::Real = 1.0; aspect::Real = 0.55)
    w = PAPER_LINEWIDTH_UNITS * fraction
    return (round(Int, w), round(Int, w * aspect))
end

"""
    paper_theme(; fontsize = 12) -> Theme

Shared manuscript theme. `fontsize` 12 Makie units = 9 pt on the page (ticks
one step smaller). Colormap `:viridis` everywhere; palette `PAPER_PALETTE`.
"""
function paper_theme(; fontsize::Real = 12)
    small = round(Int, fontsize * 0.85)   # ticks/legends: ~8 pt on the page
    Theme(
        fontsize = fontsize,
        palette = (color = PAPER_PALETTE,),
        colormap = :viridis,
        Axis = (
            titlesize = fontsize,
            titlefont = :bold,
            xlabelsize = fontsize,
            ylabelsize = fontsize,
            xticklabelsize = small,
            yticklabelsize = small,
            xgridvisible = true,
            ygridvisible = true,
            xgridcolor = (:black, 0.08),
            ygridcolor = (:black, 0.08),
            spinewidth = 0.8,
            xtickwidth = 0.8,
            ytickwidth = 0.8,
        ),
        Axis3 = (
            titlesize = fontsize,
            xlabelsize = small,
            ylabelsize = small,
            zlabelsize = small,
            xticklabelsize = small,
            yticklabelsize = small,
            zticklabelsize = small,
        ),
        Legend = (
            labelsize = small,
            titlesize = small,
            framevisible = false,
            padding = (4, 4, 2, 2),
        ),
        Colorbar = (labelsize = small, ticklabelsize = small, size = 8),
        Lines = (linewidth = 1.6,),
        Scatter = (markersize = 6,),
    )
end
