"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

include("common.jl")

# Load whole-GBR benthic and geomorphic data
gbr_benthic_path = "$(MPA_DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
gbr_benthic = Raster(gbr_benthic_path, crs=EPSG(4326), lazy=true)


gbr_morphic_path = "$(MPA_DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
gbr_geomorphic = Raster(gbr_morphic_path, crs=EPSG(4326), lazy=true)

# Get polygons of management areas - crs:GDA 2020
region_path = joinpath(MPA_DATA_DIR, "../GDA-2020/Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg")
regions_GDA2020 = GDF.read(region_path)
rename!(regions_GDA2020, Dict(:SHAPE => :geometry))

# Reproject management areas to EPSG:4326 for masking of benthic and geomorphic rasters
regions_4326 = GDF.read(region_path)
rename!(regions_4326, Dict(:SHAPE => :geometry))
regions_4326.geometry = AG.reproject(regions_4326.geometry, crs(regions_4326[1,:geometry]), crs(gbr_benthic), order=:trad)

# Reproject all raster files to GDA 2020
# Avoids the need to repeatedly reproject within each step
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_gda2020 = occursin.(reg[1:3], regions_GDA2020.AREA_DESCR)
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    # Only create files that are needed
    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"))
        target_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        target_bathy = Raster(target_bathy_path, mappedcrs=EPSG(4326), lazy=true)
        target_bathy = Rasters.resample(target_bathy; crs=crs(regions_GDA2020[1,:geometry]))

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"))
        target_slope_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        target_slope = Raster(target_slope_path, mappedcrs=EPSG(4326), lazy=true)
        target_slope = Rasters.resample(target_slope; crs=crs(regions_GDA2020[1,:geometry]))

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"), target_slope; force=true)

        target_slope = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"))
        # Trim raster to region of interest
        target_benthic = Rasters.trim(mask(gbr_benthic; with=regions_4326[reg_idx_4326, :]))
        # Reproject raster to GDA 2020
        target_benthic = resample(target_benthic; crs=crs(regions_GDA2020[1,:geometry]))

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"), target_benthic; force=true)

        target_benthic = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"))
        # Trim raster to region of interest
        target_geomorphic = Rasters.trim(mask(gbr_geomorphic; with=regions_4326[reg_idx_4326, :]))
        # Reproject raster to GDA 2020
        target_geomorphic = resample(target_geomorphic; crs=crs(regions_GDA2020[1,:geometry]))

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"), target_geomorphic; force=true)

        target_geomorphic = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"))
        # Load original bathy data to be template for wave data
        src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

        target_waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))
        target_waves_Hs = Raster(target_waves_Hs_path, key=:Hs90, crs=crs(src_bathy), mappedcrs=EPSG(4326), lazy=true)

        # Extend bounds of wave data to match bathymetry if needed
        if size(src_bathy) !== size(target_waves_Hs)
            target_waves_Hs = extend(crop(target_waves_Hs; to=src_bathy); to=AG.extent(src_bathy))
            @assert size(src_bathy) == size(target_waves_Hs)

            replace_missing!(target_waves_Hs, -9999.0)
        end

        tmp_Hs = copy(src_bathy)

        # Replace data (important: flip the y-axis!)
        tmp_Hs.data .= coalesce.(target_waves_Hs.data[:, end:-1:1], -9999.0)
        target_waves_Hs = tmp_Hs

        # Set to known missing value
        target_waves_Hs.data[target_waves_Hs.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves_Hs, -9999.0)

        # Reproject raster to GDA 2020
        target_waves_Hs = Rasters.resample(target_waves_Hs; crs=crs(regions_GDA2020[1,:geometry]))

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"), target_waves_Hs; force=true)

        src_bathy = nothing
        target_waves_Hs = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"))
        # Load original bathy data to be template for wave data
        src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

        target_waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Tp", reg)))
        target_waves_Tp = Raster(target_waves_Tp_path, key=:Tp90, crs=crs(src_bathy), mappedcrs=EPSG(4326), lazy=true)

        # Extend bounds of wave data to match bathymetry if needed
        if size(src_bathy) !== size(target_waves_Tp)
            target_waves_Tp = extend(crop(target_waves_Tp; to=src_bathy); to=AG.extent(src_bathy))
            @assert size(src_bathy) == size(target_waves_Tp)

            replace_missing!(target_waves_Tp, -9999.0)
        end

        tmp_Tp = copy(src_bathy)

        # Replace data (important: flip the y-axis!)
        tmp_Tp.data .= coalesce.(target_waves_Tp.data[:, end:-1:1], -9999.0)
        target_waves_Tp = tmp_Tp

        # Set to known missing value
        target_waves_Tp.data[target_waves_Tp.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves_Tp, -9999.0)

        # Reproject raster to GDA 2020
        target_waves_Tp = Rasters.resample(target_waves_Tp; crs=crs(regions_GDA2020[1,:geometry]))

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"), target_waves_Tp; force=true)

        src_bathy = nothing
        target_waves_Tp = nothing
        GC.gc()
    end
end
