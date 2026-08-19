"""
    Benchmark Outcome Plots (bead qjf0)

Visualization for the experiment-quality metrics that previously lived only as
JSON/CSV/Markdown tables:

- `plot_verdict_shares`: per_axis / global / tie winner shares from the
  per-axis cut counterfactual campaign (vh7e), stratified across arbitrary
  config columns (family × base_degree × budget × …), one panel per family.
- `plot_predicate_pareto`: mean_evals vs mean_recovery scatter for predicate
  threshold sweeps (ehaj.6), Pareto front highlighted, colored by
  pareto_score.

Both functions consume a `DataFrame` the caller assembled from the aggregator
outputs (verdicts.jsonl rows / analyze_predicate_sweep.jl summary.csv) —
globtimplots does no file parsing and no statistics; see
experiments/sandbox/plot_benchmark_outcomes.jl for the loading side.
"""

using Printf

const _VERDICT_ORDER = ["per_axis", "global", "tie"]
const _VERDICT_COLORS =
    Dict("per_axis" => :dodgerblue3, "global" => :darkorange2, "tie" => :gray60)

"""
    plot_verdict_shares(df::DataFrame;
                        winner_col=:winner, family_col=:family,
                        strat_cols=[:base_degree, :budget],
                        title="Per-axis cut counterfactual — winner shares")
        -> Figure

Stacked-bar shares of counterfactual winners per stratum, one panel per value
of `family_col`. Each row of `df` is one verdict record; `winner_col` values
are `"per_axis"` / `"global"` / `"tie"`. Strata are the distinct combinations
of `strat_cols` within a family, labeled `"col=val col=val"` and annotated
with the record count `n`.
"""
function plot_verdict_shares(
    df::DataFrame;
    winner_col::Symbol = :winner,
    family_col::Symbol = :family,
    strat_cols::Vector{Symbol} = [:base_degree, :budget],
    title::String = "Per-axis cut counterfactual — winner shares",
)
    nrow(df) > 0 || error("plot_verdict_shares: empty DataFrame")
    families = sort(unique(string.(df[!, family_col])))

    fig = Figure(size = (max(700, 420 * length(families)), 520))
    Label(fig[0, 1:length(families)], title; fontsize = 18, font = :bold)

    for (fi, fam) in enumerate(families)
        sub = df[string.(df[!, family_col]).==fam, :]
        strata = sort(unique([Tuple(row[c] for c in strat_cols) for row in eachrow(sub)]))
        labels = [join(["$(c)=$(v)" for (c, v) in zip(strat_cols, s)], " ") for s in strata]

        ax = Axis(
            fig[1, fi];
            title = fam,
            xticks = (1:length(strata), labels),
            xticklabelrotation = π / 4,
            ylabel = fi == 1 ? "share of records" : "",
            limits = (nothing, (0.0, 1.12)),
        )

        for (si, s) in enumerate(strata)
            rows = sub[[Tuple(row[c] for c in strat_cols) == s for row in eachrow(sub)], :]
            n = nrow(rows)
            y0 = 0.0
            for w in _VERDICT_ORDER
                share = count(==(w), string.(rows[!, winner_col])) / n
                share == 0.0 && continue
                barplot!(
                    ax,
                    [si],
                    [share];
                    offset = [y0],
                    color = _VERDICT_COLORS[w],
                    width = 0.7,
                )
                y0 += share
            end
            text!(ax, si, 1.03; text = "n=$n", align = (:center, :bottom), fontsize = 11)
        end
    end

    Legend(
        fig[2, 1:length(families)],
        [PolyElement(color = _VERDICT_COLORS[w]) for w in _VERDICT_ORDER],
        _VERDICT_ORDER;
        orientation = :horizontal,
        framevisible = false,
    )
    return fig
end

"""
    plot_predicate_pareto(df::DataFrame;
                          x_col=:mean_evals, y_col=:mean_recovery,
                          score_col=:pareto_score, label_col=nothing,
                          front_col=nothing,
                          title="Predicate sweep — evals/recovery Pareto")
        -> Figure

Scatter of predicate-sweep trials in the (cost, error) plane — both axes
minimized — colored by `score_col`, with the Pareto front drawn as a staircase
through the non-dominated points. If `front_col` names a Bool column, the
front membership computed by the analyzer is used verbatim; otherwise the
non-dominated set is derived here purely for display. `label_col` (optional)
annotates front points.
"""
function plot_predicate_pareto(
    df::DataFrame;
    x_col::Symbol = :mean_evals,
    y_col::Symbol = :mean_recovery,
    score_col::Union{Symbol,Nothing} = :pareto_score,
    label_col::Union{Symbol,Nothing} = nothing,
    front_col::Union{Symbol,Nothing} = nothing,
    title::String = "Predicate sweep — evals/recovery Pareto",
)
    nrow(df) > 0 || error("plot_predicate_pareto: empty DataFrame")
    xs = Float64.(df[!, x_col])
    ys = Float64.(df[!, y_col])

    on_front = if front_col !== nothing
        Bool.(df[!, front_col])
    else
        # Display-level dominance: point i is on the front iff no other point
        # is <= in both coordinates and < in at least one.
        [
            !any(
                (xs[j] <= xs[i]) & (ys[j] <= ys[i]) & ((xs[j] < xs[i]) | (ys[j] < ys[i])) for j in eachindex(xs)
            ) for i in eachindex(xs)
        ]
    end

    fig = Figure(size = (760, 560))
    ax = Axis(fig[1, 1]; title = title, xlabel = string(x_col), ylabel = string(y_col))

    if score_col !== nothing
        sc = scatter!(
            ax,
            xs,
            ys;
            color = Float64.(df[!, score_col]),
            colormap = :viridis,
            markersize = 12,
        )
        Colorbar(fig[1, 2], sc; label = string(score_col))
    else
        scatter!(ax, xs, ys; color = :gray40, markersize = 12)
    end

    fidx = sort(findall(on_front); by = i -> xs[i])
    if !isempty(fidx)
        stairs!(ax, xs[fidx], ys[fidx]; step = :post, color = :crimson, linewidth = 2)
        scatter!(
            ax,
            xs[fidx],
            ys[fidx];
            color = :transparent,
            strokecolor = :crimson,
            strokewidth = 2.5,
            markersize = 16,
        )
        if label_col !== nothing
            for i in fidx
                text!(
                    ax,
                    xs[i],
                    ys[i];
                    text = string(df[i, label_col]),
                    align = (:left, :bottom),
                    offset = (6, 4),
                    fontsize = 10,
                )
            end
        end
    end
    return fig
end
