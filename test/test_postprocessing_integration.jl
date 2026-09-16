"""
    test_postprocessing_integration.jl

Test that GlobtimPlots can properly import and use types from GlobtimPostProcessing.

This verifies the correct dependency direction:
    GlobtimPostProcessing (data) ← GlobtimPlots (visualization)
"""

using Test
using Dates

@testset "GlobtimPlots ← GlobtimPostProcessing Integration" begin
    @testset "GlobtimPlots can import GlobtimPostProcessing" begin
        # Load both packages
        using GlobtimPostProcessing
        using GlobtimPlots

        # Verify GlobtimPlots loaded successfully
        @test pkgversion(GlobtimPlots) isa VersionNumber
    end

    @testset "GlobtimPlots plotting functions work with GlobtimPostProcessing types" begin
        using GlobtimPostProcessing
        using GlobtimPlots
        using DataFrames

        # Create a minimal test ExperimentResult
        test_cp_df =
            DataFrame(x1 = [1.0, 2.0], x2 = [1.5, 2.5], z = [0.1, 0.2], degree = [3, 3])

        exp_result = ExperimentResult(
            "test_experiment_001",
            Dict("test_param" => 1.0, "domain_range" => 0.5),
            ["l2_discrete", "critical_points"],
            ["l2_discrete", "critical_points", "eigenvalues"],
            test_cp_df,
            Dict("runtime_seconds" => 10.5),
            nothing,
            joinpath(tempdir(), "test_exp"),
        )

        # Verify the type works
        @test exp_result isa ExperimentResult
        @test exp_result.experiment_id == "test_experiment_001"

        # Create a test campaign
        campaign = CampaignResults(
            "test_campaign_001",
            [exp_result],
            Dict("campaign_param" => "test"),
            now(),
        )

        @test campaign isa CampaignResults
        @test length(campaign.experiments) == 1

        # The campaign plots accept the real ExperimentResult / CampaignResults types
        stats = Dict(
            "approximation_quality" =>
                Dict("degrees" => [3, 4, 5], "l2_errors" => [0.3, 0.1, 0.03]),
        )
        @test create_experiment_plots(exp_result, stats; backend = Static) isa Figure
        # no params_dict in the metadata: labels fall back to positional exp_<i>
        @test generate_experiment_labels(campaign) == ["exp_1"]
        campaign_stats = Dict("test_experiment_001" => stats)
        @test create_campaign_comparison_plot(
            campaign,
            campaign_stats;
            backend = Static,
        ) isa Figure
    end
end
