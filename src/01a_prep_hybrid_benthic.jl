"""
Creating a hybrid data product combining UQ-GBRMPA dataset with ACA to gap fill missing
(in-shore) areas of the GBR.

This script creates the hybridization of the two data products for benthic data.
"""

include("common.jl")

function prep_hybrid_benthic()::Nothing
    return prep_hybrid_layer(
        "benthic",
        joinpath("benthic", "GBR10 GBRMP Benthic.tif"),
        "Benthic-Map",
        MPA_FULL_BENTHIC_IDS
    )
end
