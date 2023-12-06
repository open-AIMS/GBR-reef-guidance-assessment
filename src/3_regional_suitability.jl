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

region_features.n_slope_components .= 0
region_features.n_potential_slope .= 0


@showprogress dt=10 for reg in REGIONS
    # Partial match on region descriptor
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    d = Raster(
        joinpath(RESULT_DIR, "$(reg)_grouped_flats_95.tif"), 
        lazy=true
    )
    suitable_flats = countmap(read(d))  # get number of pixels per cluster
    region_features[reg_idx, :n_flat_components] .= length(keys(suitable_flats))
    region_features[reg_idx, :n_potential_flat] .= sum(values(suitable_flats))

    d = Raster(
        joinpath(RESULT_DIR, "$(reg)_grouped_slopes_95.tif"),
        lazy=true
    )
    suitable_slopes = countmap(read(d))  # get number of pixels per cluster
    region_features[reg_idx, :n_slope_components] .= length(keys(suitable_slopes))
    region_features[reg_idx, :n_potential_slope] .= sum(values(suitable_slopes))
end


# There's some issue with the FID column preventing successful write out
# so we ignore the column when writing results out.
GDF.write(
    joinpath(QGIS_DIR, "regional_suitability.gpkg"), 
    region_features[:, Not(:FID)];
    layer_name="region_suitability",
    geom_columns=(:geometry,)
)

subdf = region_features[:, [:AREA_DESCR, :n_flat_components, :n_potential_flat, :n_slope_components, :n_potential_slope]]
CSV.write("../qgis/potential_areas.csv", subdf)
