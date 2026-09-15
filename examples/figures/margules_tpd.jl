# margules_tpd.jl
# Phase-stability showcase figure (beads orxx.7.7 / 3ztw): the Gibbs tangent-plane
# distance of a two-suffix Margules binary, with every stationary point recovered by
# the Globtim pipeline and checked against the exact set.
#
# Run with:
#   julia --project=profiles/viz pkg/globtimplots/examples/figures/margules_tpd.jl
#   julia --project=profiles/viz pkg/globtimplots/examples/figures/margules_tpd.jl --selftest
#
# MargulesTPD is Globtim's Clapeyron-free stand-in for the TPD family, so this runs on
# profiles/viz with no thermodynamic package in the stack. The exact stationary-point
# set lives in Globtim.MARGULES_TPD_KNOWN_CPS (BigFloat Newton, 16 significant digits),
# which is what makes this a *check* and not just a picture: --selftest fails if the
# pipeline loses a stationary point or if the Poincaré–Hopf count breaks.
#
# Output: docs/showcase/assets/phase_stability_margules.png

using Globtim
using CairoMakie
using DynamicPolynomials: @polyvar
using HomotopyContinuation
using ForwardDiff
using Printf

const OUT = abspath(joinpath(@__DIR__, "..", "..", "..", "..",
                             "docs", "showcase", "assets", "phase_stability_margules.png"))
const DEGREES = [6, 8, 10, 12]
# Raw solutions of ∇p = 0 are critical points of the APPROXIMANT, displaced from the
# true ones by roughly the approximation error — measured at 3e-2 here at degree 12.
# The raw radius must therefore be loose; refinement is what makes the match tight.
const MATCH_TOL = 0.15        # logit radius for matching a RAW pipeline CP to truth
const REFINED_TOL = 1e-9      # after Newton on the true gradient

f(w) = MargulesTPD(w)
dtpd(w) = ForwardDiff.derivative(t -> MargulesTPD([t]), w)
d2tpd(w) = ForwardDiff.derivative(dtpd, w)

"""
    recover(degree) -> Vector{Float64}

Run the pipeline once at `degree` on [-5, 5] and return the recovered stationary points
in the logit coordinate, sorted.
"""
function recover(degree::Int)
    TR = TestInput(f; dim = 1, center = [0.0], GN = 400, sample_range = 5.0)
    pol = Constructor(TR, degree, basis = :chebyshev, precision = RationalPrecision)
    @polyvar x[1:1]
    real_pts = solve_polynomial_system(
        x, 1, degree, pol.coeffs;
        basis = pol.basis, precision = pol.precision, normalized = pol.normalized,
    )
    df = process_crit_pts(real_pts, f, TR)
    ws = sort(Float64.(df.x1))
    return ws, pol.nrm
end

"""
    classify(w) -> Symbol

Second-derivative classification of a stationary point of the true objective.
"""
classify(w) = d2tpd(w) > 0 ? :minimum : :maximum

"""
    refine(w0; iters = 60) -> Float64

Newton on the TRUE gradient ∂tpd/∂w, started from a critical point of the
approximant. This is the refinement step of the pipeline, done exactly here because
the objective is closed-form and 1-D.
"""
function refine(w0::Float64; iters::Int = 60)
    w = w0
    for _ in 1:iters
        g, h = dtpd(w), d2tpd(w)
        abs(h) < 1e-14 && break
        step = g / h
        w -= step
        abs(step) < 1e-15 && break
    end
    return w
end

function main()
    selftest = "--selftest" in ARGS
    known = Globtim.MARGULES_TPD_KNOWN_CPS
    @printf("exact stationary points (%d): %s\n", length(known),
            join((@sprintf("w=%.6f/%s", c.w, c.kind) for c in known), "  "))

    results = Tuple{Int,Vector{Float64},Float64,Int}[]
    for d in DEGREES
        ws, nrm = recover(d)
        matched = count(c -> any(w -> abs(w - c.w) < MATCH_TOL, ws), known)
        push!(results, (d, ws, nrm, matched))
        @printf("  degree %2d: %d raw CPs, %d/%d exact CPs matched, L2 fit err %.3e\n",
                d, length(ws), matched, length(known), nrm)
    end

    # Best arm: the lowest degree that captures the full set.
    idx = findfirst(r -> r[4] == length(known), results)
    idx === nothing && error("no degree in $DEGREES recovered all $(length(known)) stationary points")
    deg, ws, nrm, _ = results[idx]
    @printf("full capture at degree %d (L2 fit err %.3e)\n", deg, nrm)

    # Poincaré–Hopf on the recovered set, matched back to truth.
    kinds = [classify(c.w) for c in known]
    nmin = count(==(:minimum), kinds)
    nmax = count(==(:maximum), kinds)
    @printf("Poincaré–Hopf: N_min - N_max = %d - %d = %d (expected 1)\n", nmin, nmax, nmin - nmax)
    nmin - nmax == 1 || error("Poincaré–Hopf identity violated: $(nmin - nmax) != 1")

    refined = Float64[]
    raw_err, ref_err = 0.0, 0.0
    for c in known
        near = filter(w -> abs(w - c.w) < MATCH_TOL, ws)
        isempty(near) && error("stationary point w=$(c.w) not recovered at degree $deg")
        w_raw = near[1]
        w_ref = refine(w_raw)
        push!(refined, w_ref)
        raw_err = max(raw_err, abs(w_raw - c.w))
        ref_err = max(ref_err, abs(w_ref - c.w))
        @printf("  w*=%+.12f  raw %+.12f (|Δ|=%.2e)  refined %+.12f (|Δ|=%.2e)  %s\n",
                c.w, w_raw, abs(w_raw - c.w), w_ref, abs(w_ref - c.w), c.kind)
    end
    @printf("max displacement: raw %.2e, refined %.2e\n", raw_err, ref_err)
    ref_err < REFINED_TOL ||
        error("refinement did not converge to the exact set: max |Δ| = $ref_err")

    if selftest
        println("selftest OK — full stationary set recovered and Poincaré–Hopf holds")
        return
    end
    render(deg, ws, nrm, known, raw_err, ref_err)
end

function render(deg, ws, nrm, known, raw_err, ref_err)
    fig = Figure(size = (900, 560), fontsize = 15)
    ax = Axis(fig[1, 1];
        xlabel = "logit coordinate  w = ln(y₁/y₂)",
        ylabel = "tangent-plane distance  tpd(w)",
        title = "Margules binary (A = 3, z = (0.4, 0.6)) — every TPD stationary point, degree $deg")

    wgrid = range(-5, 5; length = 800)
    lines!(ax, wgrid, [f([w]) for w in wgrid]; color = :steelblue, linewidth = 2.5)
    hlines!(ax, [0.0]; color = :gray50, linestyle = :dash, linewidth = 1)

    # recovered first (large hollow), truth on top (small filled) — overlap is the point
    scatter!(ax, ws, [f([w]) for w in ws];
        color = :white, strokecolor = :crimson, strokewidth = 2.5, markersize = 20)
    for c in known
        scatter!(ax, [c.w], [c.tpd];
            color = c.kind == :minimum ? :seagreen : :darkorange,
            marker = c.kind == :minimum ? :circle : :utriangle, markersize = 11)
    end

    text!(ax, -4.7, -0.19;
        text = "N_min − N_max = 2 − 1 = 1\nL² fit error $(@sprintf("%.1e", nrm))\ndisplacement: raw $(@sprintf("%.0e", raw_err)) → refined $(@sprintf("%.0e", ref_err))",
        align = (:left, :bottom), fontsize = 13)

    elems = [
        LineElement(color = :steelblue, linewidth = 2.5),
        MarkerElement(color = :white, strokecolor = :crimson, strokewidth = 2.5,
                      marker = :circle, markersize = 14),
        MarkerElement(color = :seagreen, marker = :circle, markersize = 11),
        MarkerElement(color = :darkorange, marker = :utriangle, markersize = 11),
    ]
    labels = ["tpd(w)", "recovered by Globtim", "exact minimum", "exact maximum (feed)"]
    Legend(fig[2, 1], elems, labels; orientation = :horizontal, framevisible = false)

    mkpath(dirname(OUT))
    CairoMakie.save(OUT, fig; px_per_unit = 2)
    @printf("wrote %s\n", OUT)
end

main()
