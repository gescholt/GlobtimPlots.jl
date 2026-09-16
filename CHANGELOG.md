# Changelog

All notable changes to GlobtimPlots.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - unreleased

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
