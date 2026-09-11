# Changelog

All notable changes to GlobtimPlots.jl are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-09-10

### Added

- **Benchmark-outcome plots** — `plot_verdict_shares` and
  `plot_predicate_pareto`, for counterfactual verdict shares and the
  predicate Pareto front over threshold sweeps.
- **`paper_theme` / `paper_figsize`** — author figures at true page size so
  they need no rescaling when included, with fonts matching a Latin Modern
  manuscript body.

### Changed

- **Relicensed from GPL-3.0 to MIT.** The published `LICENSE` was already MIT;
  the README said GPL-3.0, and now agrees with it.
- Comments and docstrings no longer cite internal issue ids or monorepo-only
  script paths.

## [0.1.1] - 2026-04-30

Initial public release.

### Added

- Visualization layer for [Globtim.jl](https://github.com/gescholt/Globtim.jl) and [GlobtimPostProcessing.jl](https://github.com/gescholt/GlobtimPostProcessing.jl), built on the Makie ecosystem.
- **Critical-point and refinement plots** — scatter, level-set overlays, gradient-norm validation visualizations.
- **Convergence diagnostics** — L2 error vs. polynomial degree, residual decomposition, sample-reuse provenance.
- **Campaign comparisons** — multi-run side-by-side rendering for parameter-sweep studies.
- **Domain visualizations** — subdivision tree rendering, anisotropic-grid overlays, refinement boxes.
- **Backend extensions** (weakdeps):
  - `GlobtimDataExt` (`CSV`) — DataFrame loaders for campaign result imports.
  - `GlobtimGLMakieExt` (`GLMakie` + `DynamicPolynomials` + `Parameters`) — interactive 3D plots and polynomial visualization.
  - `GlobtimWGLMakieExt` (`WGLMakie`) — browser-rendered interactive plots for share-able UIs.

### Notes

- `CairoMakie` is the default backend (publication-quality static figures).
- `GLMakie` and `WGLMakie` are weak dependencies. Load them to activate interactive backends.
- Julia 1.12+ required.
