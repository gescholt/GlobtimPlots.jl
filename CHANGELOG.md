# Changelog

All notable changes to GlobtimPlots.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-11-19

### Changed (breaking)

- **No Makie backend is bundled.** GlobtimPlots depends on `Makie` instead of
  `CairoMakie` and re-exports Makie (`Figure`, `Axis`, `save`, ...). Load a backend
  to display or save figures: `using GlobtimPlots, CairoMakie`. Without one, `save`
  stops with Makie's "No backend available!" error. `save_plot`'s `:pdf` and
  `:svg` formats need CairoMakie. Function names are unchanged
  (`cairo_plot_polyapprox_levelset` keeps its name).
- **`plot_subdivision_tree` needs GraphMakie.** It moved to the package extension
  `GlobtimPlotsGraphMakieExt` (triggers GraphMakie, Graphs, NetworkLayout; the last
  two load with GraphMakie): add `using GraphMakie`. `TreeVizStyle` and
  `print_tree_summary` stay in the core.
- **Capture and trajectory plots need GlobtimPostProcessing.**
  `plot_capture_convergence`, `plot_capture_sparsification_combined` and
  `plot_refinement_trajectories!` moved to the package extension
  `GlobtimPlotsPostProcessingExt`: add `using GlobtimPostProcessing`.
  `RefinementTrajectoryStyle` stays in the core. GlobtimPostProcessing is no longer
  a dependency.
- Without its package, each moved function raises an error naming the package to
  load ("... requires GraphMakie.jl" / "... requires GlobtimPostProcessing.jl").

### Removed

- Dependencies `CairoMakie`, `GraphMakie`, `Graphs`, `NetworkLayout`,
  `GlobtimPostProcessing` (now weak dependencies), `GeometryBasics` and
  `LaTeXStrings` (`Rect`, `Point2f` and `LaTeXString` come from Makie).

### Fixed

- Compat for GlobtimPostProcessing admits 0.2.

## [0.1.0] - 2026-11-19

Initial public release: the visualization layer for
[Globtim.jl](https://github.com/gescholt/Globtim.jl) and
[GlobtimPostProcessing.jl](https://github.com/gescholt/GlobtimPostProcessing.jl),
built on the Makie ecosystem (CairoMakie by default, GLMakie / WGLMakie via
package extensions).

### Added

- **Level sets** — `cairo_plot_polyapprox_levelset` (approximant contours with
  raw / refined critical points), `plot_polyapprox_3d` (surface, GLMakie or
  WGLMakie extension), 3D level-set slider and animation helpers.
- **Morse diagnostics** — Hessian eigenvalue, condition-number and norm plots per
  critical point, colored by classification.
- **Degree sweeps** — discrete L2 error, refinement-displacement and captured-point
  statistics against the polynomial degree.
- **Adaptive subdivision** — partition, objective overlay, per-cell error heatmap,
  tree diagram, and the GLMakie zoom-refine / click-to-optimize explorers.
- **Experiments and campaigns** — label-aware summaries of tracked statistics for
  one experiment or a campaign (`create_experiment_plots`,
  `create_campaign_comparison_plot`).
- **Capture analysis, refinement trajectories, benchmark-outcome plots**
  (`plot_capture_convergence`, `plot_refinement_trajectories!`,
  `plot_verdict_shares`, `plot_predicate_pareto`).
- **Backend extensions** (weakdeps): `GlobtimPlotsGLMakieExt` (`GLMakie`) and
  `GlobtimPlotsWGLMakieExt` (`WGLMakie`).

### Notes

- Julia 1.12+ required. Licensed under MIT.
- Research-specific figure code that lived in this package while it was
  developed inside a monorepo (RL training dashboards, Lotka-Volterra 4D
  experiment plots, a manuscript figure theme) is not part of the release.
