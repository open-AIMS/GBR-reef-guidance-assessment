"""
Prepare data for analysis.

Crop GBR-wide GBRMPA rasters into management regions.
Reproject all data from WGS84 / UTM Zone 54 - 56 into consistent crs GDA-2020.
Ensure all rasters are the same size/shape for each region of interest.
"""

include("common.jl")

# 1. Processing of GBRMPA zoning geopackage to select only zones for site exclusion
if !isfile(joinpath(MPA_OUTPUT_DIR, "GBRMPA_zone_exclusion.gpkg"))
    GBRMPA_zoning_poly = GDF.read(joinpath(GDA2020_DATA_DIR, "Great_Barrier_Reef_Marine_Park_Zoning_20_4418126048110066699.gpkg"))
    GBRMPA_zoning_poly = GBRMPA_zoning_poly[GBRMPA_zoning_poly.TYPE .∈ [MPA_EXCLUSION_ZONES], :]
    rename!(GBRMPA_zoning_poly, :SHAPE => :geometry)

    GDF.write(joinpath(MPA_OUTPUT_DIR, "GBRMPA_zone_exclusion.gpkg"), GBRMPA_zoning_poly; crs=EPSG(7844))
end

# 2. Processing of MPA files into smaller GDA-2020 rasters
# Loading regions_4326 for cropping of vector and raster data.
regions_4326 = GDF.read(REGION_PATH_4326)

# Loading GBR-wide data
gbr_benthic_path = "$(MPA_DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
gbr_benthic = Raster(gbr_benthic_path, crs=EPSG(4326), lazy=true)

gbr_morphic_path = "$(MPA_DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
gbr_geomorphic = Raster(gbr_morphic_path, crs=EPSG(4326), lazy=true)

aca_turbid_path = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
aca_turbid = Raster(aca_turbid_path, mappedcrs=EPSG(4326), lazy=true)

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"))
        target_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        target_bathy = Raster(target_bathy_path, mappedcrs=EPSG(4326), lazy=true)
        target_bathy = Rasters.resample(target_bathy; crs=GDA2020_crs)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"))
        target_slope_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        target_slope = Raster(target_slope_path, mappedcrs=EPSG(4326), lazy=true)
        target_slope = Rasters.resample(target_slope; crs=GDA2020_crs)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"), target_slope; force=true)

        target_slope = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"))
        bathy_gda2020_path = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif")
        bathy_gda2020 = Raster(bathy_gda2020_path; crs=EPSG(7844), lazy=true)

        # Trim raster to region of interest
        target_benthic = Rasters.trim(mask(gbr_benthic; with=regions_4326[reg_idx_4326, :]))
        # Reproject raster to GDA 2020
        target_benthic = resample(target_benthic; to=bathy_gda2020)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"), target_benthic; force=true)

        bathy_gda2020 = nothing
        target_benthic = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"))
        bathy_gda2020_path = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif")
        bathy_gda2020 = Raster(bathy_gda2020_path; crs=EPSG(7844), lazy=true)

        # Trim raster to region of interest
        target_geomorphic = Rasters.trim(mask(gbr_geomorphic; with=regions_4326[reg_idx_4326, :]))
        # Reproject raster to GDA 2020
        target_geomorphic = resample(target_geomorphic; to=bathy_gda2020)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"), target_geomorphic; force=true)

        bathy_gda2020 = nothing
        target_geomorphic = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"))
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
        tmp_Hs.data .= target_waves_Hs.data[:, end:-1:1]
        target_waves_Hs = tmp_Hs

        # Set to known missing value
        target_waves_Hs.data[target_waves_Hs.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves_Hs, -9999.0)

        # Reproject raster to GDA 2020
        target_waves_Hs = Rasters.resample(target_waves_Hs; crs=GDA2020_crs)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"), target_waves_Hs; force=true)

        src_bathy = nothing
        target_waves_Hs = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"))
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
        tmp_Tp.data .= target_waves_Tp.data[:, end:-1:1]
        target_waves_Tp = tmp_Tp

        # Set to known missing value
        target_waves_Tp.data[target_waves_Tp.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves_Tp, -9999.0)

        # Reproject raster to GDA 2020
        target_waves_Tp = Rasters.resample(target_waves_Tp; crs=GDA2020_crs)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"), target_waves_Tp; force=true)

        src_bathy = nothing
        target_waves_Tp = nothing
        GC.gc()
    end

    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"))
        bathy_gda2020_path = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif")
        bathy_gda2020 = Raster(bathy_gda2020_path; crs=EPSG(7844), lazy=true)

        target_turbid = Rasters.crop(aca_turbid; to=regions_4326[reg_idx_4326, :])
        target_turbid = Rasters.trim(mask(target_turbid; with=regions_4326[reg_idx_4326, :]))
        target_turbid = resample(target_turbid; to=bathy_gda2020)

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"), target_turbid; force=true)

        bathy_gda2020 = nothing
        target_turbid = nothing
        GC.gc()
    end

    if reg == "Townsville-Whitsunday"
        if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"))
            tsv_rugosity_path = joinpath(RUG_DATA_DIR, "std25_Rugosity_Townsville-Whitsunday.tif")
            tsv_rugosity = Raster(tsv_rugosity_path; lazy=true)
            tsv_rugosity = Rasters.resample(tsv_rugosity; crs=GDA2020_crs)

            write(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"), tsv_rugosity; force=true)

            tsv_rugosity = nothing
            GC.gc()
        end
    end
end
