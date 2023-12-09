"""
Collate number of potentially suitable locations per reef, as defined by
GBRMPA Features.
"""

using Rasters
import GeoDataFrames as GDF
import ArchGDAL as AG
using CSV

using Statistics, StatsBase
using Glob
using ProgressMeter

include("common.jl")

# D:\development\ADRIA_data\spatial_datasets\Bathy data 10m\features
reef_path = joinpath(
    DATA_DIR, 
    "features", 
    "Great_Barrier_Reef_Features.shp"
)

reef_features = GDF.read(reef_path)
reef_features.n_potential_flats .= 0
reef_features.n_potential_slopes .= 0
reef_features.n_potential_flats_ha .= 0.0
reef_features.n_potential_slopes_ha .= 0.0

reef_features.flats_score .= 0.0
reef_features.slopes_score .= 0.0

"""
    count_suitable(raster, reef)::Int64

Count number of suitable pixels in an area as defined by the `reef` geometry, for the
given `raster`.

This has to be applied one by one for each reef as any malformed geometries lead to a crash.
"""
function count_suitable(raster, reef)::Int64
    local total::Union{Int, Missing} = 0
    try
        total = Rasters.zonal(sum, raster; of=reef, shape=:polygon, boundary=:touches)
        if ismissing(total)
            total = 0
        end

    catch err
        if !(err isa TaskFailedException)
            rethrow(err)
        end

        @info "Failed to extract stats for $(reef.GBR_NAME) - $(reef.UNIQUE_ID)"
        total = 0
    end

    return total
end

# Loop over each reef and count number of slopes and flats that meet criteria
@showprogress dt=10 desc="Collating zonal stats..." for reg in REGIONS
    # Load raster
    target_flats = Raster(
        joinpath(RESULT_DIR, "$(reg)_suitable_flats_85.tif"),
        mappedcrs=EPSG(4326)
    )

    target_slopes = Raster(
        joinpath(RESULT_DIR, "$(reg)_suitable_slopes_85.tif"),
        mappedcrs=EPSG(4326)
    )

    # For some reason the geometries references in the GeoDataFrame are referenced
    # so any transforms also affect the original.
    # Any later transforms appear invalid (as it is no longer in EPSG:4326) and so leads
    # to a crash. 
    # Taking a copy does not work, so a quick workaround is to simply read the geometries
    # in again.
    reefs = GDF.read(reef_path)
    AG.reproject(reefs.geometry, GFT.EPSG(4326), crs(target_flats); order=:trad)

    for reef in eachrow(reefs)
        # Count number of locations that meet flats and slopes criteria

        # Match by OBJECTID, which should be the row number but don't trust it...
        target_row = reef_features.OBJECTID .== reef.OBJECTID
        reef_features[target_row, :n_potential_flats] .= count_suitable(target_flats, reef.geometry)
        reef_features[target_row, :n_potential_flats_ha] .= reef_features[target_row, :n_potential_flats] / 100.0

        # Do again for slopes
        reef_features[target_row, :n_potential_slopes] .= count_suitable(target_slopes, reef.geometry)
        reef_features[target_row, :n_potential_slopes_ha] .= reef_features[target_row, :n_potential_slopes] / 100.0

        reef_features[target_row, :flats_score] .= reef_features[target_row, :n_potential_flats_ha] / reef_features[target_row, :Area_HA]
        reef_features[target_row, :slopes_score] .= reef_features[target_row, :n_potential_slopes_ha] / reef_features[target_row, :Area_HA]
    end

    target_flats = nothing
    target_slopes = nothing
    GC.gc()
end

reef_features[:, :flats_score] .= reef_features[:, :flats_score] / maximum(reef_features[:, :flats_score])
reef_features[:, :slopes_score] .= reef_features[:, :slopes_score] / maximum(reef_features[:, :slopes_score])

GDF.write(
    joinpath(QGIS_DIR, "reef_suitability.shp"), 
    reef_features[:, [:geometry, :LOC_NAME_S, :UNIQUE_ID, :Area_HA, :n_potential_flats, :n_potential_flats_ha, :n_potential_slopes, :n_potential_slopes_ha, :flats_score, :slopes_score]];
    layer_name="reef_suitability",
    geom_columns=(:geometry,),
    crs=EPSG(4326)
)

subdf = reef_features[:, [:LOC_NAME_S, :UNIQUE_ID, :Area_HA, :n_potential_flats, :n_potential_flats_ha, :n_potential_slopes, :n_potential_slopes_ha, :flats_score, :slopes_score]]
CSV.write("../qgis/potential_reef_areas.csv", subdf)
