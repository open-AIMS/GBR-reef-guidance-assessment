"""Count potential locations for each region.

This is a bit suspect, and superceded by Step 4.
Kept as it may come in handy later.
"""

using Rasters
import GeoDataFrames as GDF
import ArchGDAL as AG
using CSV

using Statistics, StatsBase
using Glob

include("common.jl")

region_path = joinpath(
    DATA_DIR, 
    "zones", 
    "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
)

region_features = copy(GDF.read(region_path))
region_features.n_flat_components .= 0
region_features.n_potential_flat .= 0
region_features.n_potential_flat_ha .= 0.0

region_features.n_slope_components .= 0
region_features.n_potential_slope .= 0
region_features.n_potential_slope_ha .= 0.0


@showprogress dt=10 for reg in REGIONS
    # Partial match on region descriptor
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    flats = Raster(
        joinpath(RESULT_DIR, "$(reg)_grouped_flats_85.tif")
    )
    grouped_flats = countmap(flats)  # get number of pixels per cluster
    region_features[reg_idx, :n_flat_components] .= length(keys(grouped_flats))
    region_features[reg_idx, :n_potential_flat] .= sum(values(grouped_flats))
    region_features[reg_idx, :n_potential_flat_ha] .= count(flats .> 0.0) / 100.0
    flats = nothing
    GC.gc()

    slopes = Raster(
        joinpath(RESULT_DIR, "$(reg)_grouped_slopes_85.tif")
    )
    grouped_slopes = countmap(slopes)  # get number of pixels per cluster
    region_features[reg_idx, :n_slope_components] .= length(keys(grouped_slopes))
    region_features[reg_idx, :n_potential_slope] .= sum(values(grouped_slopes))
    region_features[reg_idx, :n_potential_slope_ha] .= count(slopes .> 0.0) / 100.0
    slopes = nothing
    GC.gc()
end

# There's some issue with the FID column preventing successful write out
# so we ignore the column when writing results out.
GDF.write(
    joinpath(QGIS_DIR, "regional_suitability.gpkg"), 
    region_features[:, Not(:FID)];
    layer_name="region_suitability",
    geom_columns=(:geometry,)
)

subdf = region_features[:, [:AREA_DESCR, :n_flat_components, :n_potential_flat_ha, :n_slope_components, :n_potential_slope_ha]]
CSV.write("../qgis/region_potential_areas.csv", subdf)
