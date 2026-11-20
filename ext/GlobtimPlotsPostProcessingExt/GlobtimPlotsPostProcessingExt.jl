module GlobtimPlotsPostProcessingExt

# GlobtimPostProcessing extension for GlobtimPlots: the plots whose inputs are
# GlobtimPostProcessing types (CaptureResult, KnownCriticalPoints,
# CriticalPointRefinementResult). Loaded by `using GlobtimPostProcessing`.

using GlobtimPlots
using GlobtimPlots: RefinementTrajectoryStyle
using GlobtimPostProcessing: CaptureResult, KnownCriticalPoints, CriticalPointRefinementResult
using Makie
using Printf
using Random: MersenneTwister
using LinearAlgebra: norm

include("capture_analysis_plots.jl")
include("refinement_trajectory_viz.jl")

end # module GlobtimPlotsPostProcessingExt
