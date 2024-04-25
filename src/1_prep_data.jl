"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

include("common.jl")


# Reproject benthic data into common CRSs
# Avoids the need to repeatedly reproject within each step

gbr_benthic_path = "$(MPA_DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
gbr_benthic = Raster(gbr_benthic_path, crs=EPSG(4326), lazy=true)


gbr_morphic_path = "$(MPA_DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
gbr_geomorphic = Raster(gbr_morphic_path, crs=EPSG(4326), lazy=true)


# Get polygon of management areas
region_path = joinpath(
    MPA_DATA_DIR,
    "zones",
    "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
)
region_features = GDF.read(region_path)

@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    if !isfile(joinpath(OUTPUT_DIR, "$(reg)_benthic.tif"))
        # Only recreate files if needed
        target_benthic = Rasters.trim(mask(gbr_benthic; with=region_features[reg_idx, :]))
        target_benthic = resample(target_benthic, to=src_bathy)
        write(joinpath(OUTPUT_DIR, "$(reg)_benthic.tif"), target_benthic; force=true)

        target_benthic = nothing
        GC.gc()
    end

    if !isfile(joinpath(OUTPUT_DIR, "$(reg)_geomorphic.tif"))
        target_geomorphic = Rasters.trim(mask(gbr_geomorphic; with=region_features[reg_idx, :]))
        target_geomorphic = resample(target_geomorphic, to=src_bathy)
        write(joinpath(OUTPUT_DIR, "$(reg)_geomorphic.tif"), target_geomorphic; force=true)

        target_geomorphic = nothing
        GC.gc()
    end

    if !isfile(joinpath(OUTPUT_DIR, "$(reg)_waves.tif"))
        target_waves_path = first(glob("*.nc", joinpath(MPA_DATA_DIR, "waves", reg)))
        target_waves = Raster(target_waves_path, key=:Hs90, crs=crs(src_bathy), mappedcrs=EPSG(4326), lazy=true)

        # Extend bounds of wave data to match bathymetry if needed
        if size(src_bathy) !== size(target_waves)
            target_waves = extend(crop(target_waves; to=src_bathy); to=AG.extent(src_bathy))
            @assert size(src_bathy) == size(target_waves)
        end

        # Hacky workaround:
        # Due to projection schenanigans, the extents appears to be offset by ~5m.
        # The extent of the data may also cross two UTM zones, even if the wave data
        # is well within the bathymetry bounds.
        # We assume the wave coordinates are incorrect and move on by copying the bathymetry
        # data structure and replace its values with wave data.
        tmp = copy(src_bathy)

        # Replace data (important: flip the y-axis!)
        tmp.data .= target_waves.data[:, end:-1:1]
        target_waves = tmp

        # Set to known missing value
        target_waves.data[target_waves.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves, -9999.0)

        write(joinpath(OUTPUT_DIR, "$(reg)_waves.tif"), target_waves; force=true)

        target_waves = nothing
        GC.gc()
    end
end
