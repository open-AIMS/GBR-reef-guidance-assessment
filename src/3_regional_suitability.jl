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
region_features.suitable_flat_area_ha .= 0
region_features.suitable_flat_area_m2 .= 0
region_features.potential_flat_area_ha .= 0
region_features.potential_flat_area_m2 .= 0

region_features.n_slope_components .= 0
region_features.suitable_slope_area_ha .= 0
region_features.suitable_slope_area_m2 .= 0
region_features.potential_slope_area_ha .= 0
region_features.potential_slope_area_m2 .= 0


@showprogress dt=10 for reg in REGIONS
    d = Raster(
        joinpath(RESULT_DIR, "$(reg)_grouped_flats_95.tif"), 
        lazy=true
    )

    # Load pre-prepared benthic data
    src_benthic_path = "../figs/$(reg)_benthic.tif"
    src_benthic = Raster(src_benthic_path, lazy=true)

    src_geomorphic_path = "../figs/$(reg)_geomorphic.tif"
    src_geomorphic = Raster(src_geomorphic_path, lazy=true)

    # Count pixels that are rock or coral/algae on flats
    potential_flat_area = count(read(
        (src_geomorphic .∈ [FLAT_IDS]) .&
        (src_benthic .∈ [BENTHIC_IDS])
    ))

    # Partial match on region descriptor
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    suitable_flats = countmap(read(d))  # get number of pixels
    region_features[reg_idx, :n_flat_components] .= length(keys(suitable_flats))
    region_features[reg_idx, :suitable_flat_area_ha] .= sum(values(suitable_flats))

    # Because each pixel is a [10m ⋅ 10m] area, multiply by 100 to get area in m² 
    region_features[reg_idx, :suitable_flat_area_m2] .= sum(values(suitable_flats)) * 100

    region_features[reg_idx, :potential_flat_area_ha] .= potential_flat_area
    region_features[reg_idx, :potential_flat_area_m2] .= potential_flat_area * 100

    d = Raster(
        joinpath(RESULT_DIR, "$(reg)_grouped_slopes_95.tif"), 
        lazy=true
    )

    suitable_slopes = countmap(read(d))  # get number of pixels
    region_features[reg_idx, :n_slope_components] .= length(keys(suitable_slopes))
    region_features[reg_idx, :suitable_slope_area_ha] .= sum(values(suitable_slopes))

    # Because each pixel is a [10m ⋅ 10m] area, multiply by 100 to get area in m² 
    region_features[reg_idx, :suitable_slope_area_m2] .= sum(values(suitable_slopes)) * 100

    potential_slope_area = count(read(
        (src_geomorphic .∈ [SLOPE_IDS]) .&
        (src_benthic .∈ [BENTHIC_IDS])
    ))

    region_features[reg_idx, :potential_slope_area_ha] .= potential_slope_area
    region_features[reg_idx, :potential_slope_area_m2] .= potential_slope_area * 100
end


# There's some issue with the FID column preventing successful write out
# so we ignore the column when writing results out.
GDF.write(
    joinpath(QGIS_DIR, "regional_suitability.gpkg"), 
    region_features[:, Not(:FID)];
    layer_name="region_suitability",
    geom_columns=(:geometry,)
)

subdf = region_features[:, [:AREA_DESCR, :n_flat_components, :potential_flat_area_ha, :suitable_flat_area_ha, :potential_slope_area_ha, :suitable_slope_area_ha]]
CSV.write("../qgis/potential_areas.csv", subdf)
