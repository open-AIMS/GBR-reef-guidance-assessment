"""
Creating a hybrid data product combining UQ-GBRMPA dataset with ACA to gap fill missing
(in-shore) areas of the GBR.

This script creates the hybridization of the two data products for geomorphic data.
"""

include("common.jl")

function prep_hybrid_geomorphic()::Nothing
    return prep_hybrid_layer(
        "geomorphic",
        joinpath("geomorphic", "GBR10 GBRMP Geomorphic.tif"),
        "Geomorphic-Map",
        MPA_GEOMORPHIC_IDS
    )
end
