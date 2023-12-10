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
reef_features.region .= ""
reef_features.reef_name .= ""
reef_features.n_flat .= 0
reef_features.n_slope .= 0
reef_features.flat_ha .= 0.0
reef_features.slope_ha .= 0.0

reef_features.flat_scr .= 0.0
reef_features.slope_scr .= 0.0

"""
    count_suitable(raster, reef)::Int64

Count number of suitable pixels in an area as defined by the `reef` geometry, for the
given `raster`.

This has to be applied one by one for each reef as any malformed geometries lead to a crash.
"""
function count_suitable(raster, reef)::Union{Int64, Missing}
    local total::Union{Int, Missing} = missing
    try
        total = Rasters.zonal(sum, raster; of=reef, shape=:polygon, boundary=:touches)
    catch err
        if !(err isa TaskFailedException)
            rethrow(err)
        end

        @info "Failed to extract stats for $(reef.GBR_NAME) - $(reef.UNIQUE_ID)"
        total = missing
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
    reefs.geometry = AG.reproject(reefs.geometry, GFT.EPSG(4326), crs(target_flats); order=:trad)

    for (target_row, reef) in enumerate(eachrow(reefs))
        # Count number of locations that meet flats and slopes criteria.
        # Have to loop over each reef individually as some geometries causes a crash.
        # These should be safe to skip, and will resolve to `missing`.
        if reef_features[target_row, :region] != ""
            # Reef already associated with a region and therefore have been parsed in a
            # previous loop, so we skip.
            continue
        end

        flat_val = count_suitable(target_flats, reef.geometry)
        if !ismissing(flat_val)
            reef_features[target_row, :n_flat] = flat_val
            reef_features[target_row, :flat_ha] = flat_val / 100.0
        end

        # Do again for slopes
        slope_val = count_suitable(target_slopes, reef.geometry)
        if !ismissing(slope_val)
            reef_features[target_row, :n_slope] = slope_val
            reef_features[target_row, :slope_ha] = slope_val / 100.0
        end

        reef_features[target_row, :reef_name] = reef.LOC_NAME_S

        if !ismissing(flat_val) || !ismissing(slope_val)
            reef_features[target_row, :region] = reg
        end
    end
    
    target_flats = nothing
    target_slopes = nothing
    GC.gc()
end

valid_locs = reef_features.Area_HA .!= 0
reef_features[valid_locs, :flat_scr] .= round.(reef_features[valid_locs, :flat_ha] ./ reef_features[valid_locs, :Area_HA], digits=4)
reef_features[valid_locs, :slope_scr] .= round.(reef_features[valid_locs, :slope_ha] ./ reef_features[valid_locs, :Area_HA], digits=4)

for reg in REGIONS
    target_reg = reef_features.region .== reg
    reef_features[target_reg, :flat_scr] .= round.(reef_features[target_reg, :flat_scr] / maximum(reef_features[target_reg, :flat_scr]), digits=4)
    reef_features[target_reg, :slope_scr] .= round.(reef_features[target_reg, :slope_scr] / maximum(reef_features[target_reg, :slope_scr]), digits=4)
end

# Have to write out results as shapefile because of ArcGIS not handling GeoPackages for
# some reason...
GDF.write(
    joinpath(QGIS_DIR, "reef_suitability.shp"),
    reef_features[:, [:geometry, :region, :reef_name, :flat_ha, :slope_ha, :Area_HA, :n_flat, :n_slope, :flat_scr, :slope_scr, :UNIQUE_ID]];
    layer_name="reef_suitability",
    geom_columns=(:geometry,),
    crs=EPSG(4326)
)

subdf = reef_features[:, [:region, :reef_name, :flat_ha, :slope_ha, :Area_HA, :n_flat, :n_slope, :flat_scr, :slope_scr, :UNIQUE_ID]]
CSV.write("../qgis/potential_reef_areas.csv", subdf)
