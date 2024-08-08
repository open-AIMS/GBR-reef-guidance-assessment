"""
Collate number of potentially suitable locations per reef, as defined by
GBRMPA Features and analyses from GBRMPA data.
"""

using CSV

include("common.jl")

# We need `Area_HA` and `UNIQUE_ID` from the GBR features dataset provided,
# however this is in crs GDA-94 so we have to reproject this dataset to GDA-2020.
reef_features = GDF.read(REEF_PATH_GDA94)
reef_features.geometry = AG.reproject(reef_features.geometry, crs(reef_features[1, :geometry]), GDA2020_crs; order=:trad)

reef_features.region .= ""
reef_features.reef_name .= ""

reef_features.n_flat .= 0
reef_features.n_slope .= 0

reef_features.flat_ha .= 0.0
reef_features.slope_ha .= 0.0

reef_features.flat_scr .= 0.0
reef_features.slope_scr .= 0.0

"""
    count_suitable(raster, reef)::Union{Int64, Missing}

Count the number of pixels in an area with greater than 95% surrounding suitability,
as defined by the `reef` geometry, for the given `raster`.

This has to be applied one by one for each reef as any malformed geometries lead to a crash.
"""
function count_suitable(raster, reef)::Union{Int64,Missing}
    local total::Union{Int,Missing} = missing
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
@showprogress dt = 10 desc = "Collating zonal stats..." for reg in REGIONS
    # Load rasters and identify cells that are greater than or equal to 95% suitability
    target_flats = Raster(
        joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_flats.tif"),
        crs=EPSG(7844),
        lazy=true
    )
    target_flats = read(target_flats .>= 95)

    target_slopes = Raster(
        joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif"),
        crs=EPSG(7844),
        lazy=true
    )
    target_slopes = read(target_slopes .>= 95)

    reefs = reef_features

    for (target_row, reef) in enumerate(eachrow(reefs))
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

# Calculate the proportion of each reef area that meets 95% suitability criteria
valid_locs = reef_features.Area_HA .!= 0
reef_features[valid_locs, :flat_scr] .= round.(reef_features[valid_locs, :flat_ha] ./ reef_features[valid_locs, :Area_HA], digits=4)
reef_features[valid_locs, :slope_scr] .= round.(reef_features[valid_locs, :slope_ha] ./ reef_features[valid_locs, :Area_HA], digits=4)

for reg in REGIONS
    # `scr` then gives an index to find the reefs within a region that have the most
    # proportion of their area suitable for deployments
    target_reg = reef_features.region .== reg
    reef_features[target_reg, :flat_scr] .= round.(reef_features[target_reg, :flat_scr] / maximum(reef_features[target_reg, :flat_scr]), digits=4)
    reef_features[target_reg, :slope_scr] .= round.(reef_features[target_reg, :slope_scr] / maximum(reef_features[target_reg, :slope_scr]), digits=4)
end

# Write data to geopackage
GDF.write(
    joinpath(MPA_QGIS_DIR, "reef_suitability.gpkg"),
    reef_features[:, [
        :geometry,
        :region,
        :reef_name,
        :flat_ha,
        :slope_ha,
        :Area_HA,
        :n_flat,
        :n_slope,
        :flat_scr,
        :slope_scr,
        :UNIQUE_ID
    ]];
    crs=EPSG(7844)
)

# Write data to csv file
subdf = reef_features[:, [
    :region,
    :reef_name,
    :flat_ha,
    :slope_ha,
    :Area_HA,
    :n_flat,
    :n_slope,
    :flat_scr,
    :slope_scr,
    :UNIQUE_ID
]]
CSV.write(joinpath(MPA_QGIS_DIR, "potential_reef_areas.csv"), subdf)

# Rank reefs by their regional suitability score
reef_scores = reef_features[:, [
    :geometry,
    :region,
    :reef_name,
    :flat_ha,
    :slope_ha,
    :Area_HA,
    :n_flat,
    :n_slope,
    :flat_scr,
    :slope_scr,
    :UNIQUE_ID
]]

# Keep reefs that have some suitability score
reefs_with_scores = reef_scores[(reef_scores.flat_scr .!== 0.0) .| (reef_scores.slope_scr .!== 0.0), :]

# Rank reefs by flat_scr and include the top 10 reefs
highest_flats = DataFrames.combine(groupby(reefs_with_scores,:region), sdf -> sort(sdf,:flat_scr; rev=true), :region => eachindex => :rank)
highest_flats[highest_flats.rank .∈ [1:10], :]
GDF.write(joinpath(MPA_QGIS_DIR, "highest_ranked_reefs_flats.gpkg"), highest_flats; crs=EPSG(7844))

# Rank reefs by slope_scr and include the top 10 reefs
highest_slopes = DataFrames.combine(groupby(reefs_with_scores, :region), sdf -> sort(sdf, :slope_scr; rev=true), :region => eachindex => :rank)
highest_slopes[highest_slopes.rank .∈ [1:10], :]
GDF.write(joinpath(MPA_QGIS_DIR, "highest_ranked_reefs_slopes.gpkg"), highest_slopes; crs=EPSG(7844))
