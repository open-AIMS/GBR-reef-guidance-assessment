"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

import GeoDataFrames as GDF
import ArchGDAL as AG
using GeoInterface

include("common.jl")

if !isfile(first(glob("*.gpkg", ACA_OUTPUT_DIR)))
    geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson"))
    benthic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Benthic-Map", "benthic.geojson"))
    reef_poly = GDF.read(joinpath(ACA_DATA_DIR, "Reef-Extent", "reefextent.geojson"))

    target_geomorphic_features = geomorphic_poly.class .∈ Ref([ACA_FLAT_IDS..., ACA_SLOPE_IDS])
    target_benthic_features = benthic_poly.class .∈ Ref(ACA_BENTHIC_IDS)

    geomorphic_poly = geomorphic_poly[target_geomorphic_features, :]
    benthic_poly = benthic_poly[target_benthic_features, :]
    GC.gc()

    # Assuming we skip any reefs where we do not have data across all criteria
    skipped_reefs = fill("", size(reef_poly, 1))
    @floop for (idx, reef) in enumerate(eachrow(reef_poly))
        reef_geom = reef.geometry
        has_geomorphic = AG.intersects.([reef_geom], geomorphic_poly.geometry)
        if all(has_geomorphic .== false)
            skipped_reefs[idx] = "No geomorphic"
            continue
        end

        has_benthic = AG.intersects.([reef_geom], benthic_poly.geometry)
        if all(has_benthic .== false)
            skipped_reefs[idx] = "No benthic"
            continue
        end
    GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_benthic_geomorphic.gpkg"), reef_poly[skipped_reefs .== "", :])end
end
