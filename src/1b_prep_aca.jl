"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

include("common.jl")

# Loading regions_4326 for cropping of vector and raster data.
regions_4326 = GDF.read(REGION_PATH_4326)

# 1. Processing of geojson files into smaller-GBRMPA regions

# Loading GBR-wide data
geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson"))
benthic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Benthic-Map", "benthic.geojson"))

target_flats = geomorphic_poly.class .∈ Ref(ACA_FLAT_IDS)
target_flat_poly = geomorphic_poly[target_flats, :]

target_slopes = geomorphic_poly.class .∈ Ref(ACA_SLOPE_IDS)
target_slope_poly = geomorphic_poly[target_slopes, :]

target_benthic_features = benthic_poly.class .∈ Ref(ACA_BENTHIC_IDS)
target_benthic_poly = benthic_poly[target_benthic_features, :]

@floop for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)
    region_4326_geom = regions_4326[reg_idx_4326, :geometry][1]

    if !isfile(joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg).gpkg"))
        flat_is_in_region = AG.contains.([region_4326_geom], target_flat_poly.geometry)
        target_flats_reg = target_flat_poly[flat_is_in_region, :]

        target_flats_reg.geometry = AG.reproject(target_flats_reg.geometry, crs(region_4326_geom), GDA2020_crs; order=:trad)

        GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg).gpkg"), target_flats_reg; crs=EPSG(7844))
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg).gpkg"))
        slope_is_in_region = AG.contains.([region_4326_geom], target_slope_poly.geometry)
        target_slopes_reg = target_slope_poly[slope_is_in_region, :]

        target_slopes_reg.geometry = AG.reproject(target_slopes_reg.geometry, crs(region_4326_geom), GDA2020_crs; order=:trad)

        GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg).gpkg"), target_slopes_reg; crs=EPSG(7844))
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg).gpkg"))
        ben_is_in_region = AG.contains.([region_4326_geom], benthic_poly.geometry)
        benthic_reg = benthic_poly[ben_is_in_region, :]

        benthic_reg.geometry = AG.reproject(benthic_reg.geometry, crs(region_4326_geom), GDA2020_crs; order=:trad)

        GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg).gpkg"), benthic_reg; crs=EPSG(7844))
    end
end

# 2. Processing of reef-wide rasters into management region scale rasters, and reprojection
# of bathymetry, turbidity and wave raster data.

# Loading GBR-wide data
aca_bathy_path = "$(ACA_DATA_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
aca_bathy = Raster(aca_bathy_path, mappedcrs=EPSG(4326), lazy=true)

aca_turbid_path = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
aca_turbid = Raster(aca_turbid_path, mappedcrs=EPSG(4326), lazy=true)

@showprogress dt = 10 "Prepping bathymetry/turbidity/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"))
        # Cropping the raster speeds up the trim(mask()) step and reduces memory use - checked output in qgis
        target_bathy = Rasters.crop(aca_bathy; to=regions_4326[reg_idx_4326, :])
        target_bathy = Rasters.trim(mask(target_bathy; with=regions_4326[reg_idx_4326, :]))
        target_bathy = Rasters.resample(target_bathy; crs=GDA2020_crs)

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"))
        # Using aca_bathy as the basis for reprojecting turbidity data to GDA2020
        bathy_gda2020_path = joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif")
        bathy_gda2020 = Raster(bathy_gda2020_path; crs=EPSG(7844), lazy=true)

        target_turbid = Rasters.crop(aca_turbid; to=regions_4326[reg_idx_4326, :])
        target_turbid = Rasters.trim(mask(target_turbid; with=regions_4326[reg_idx_4326, :]))
        target_turbid = resample(target_turbid; to=bathy_gda2020)

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"), target_turbid; force=true)

        target_turbid = nothing
        src_aca_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"))
        # Need MPA bathy data as a template for wave data due to unknown crs/format
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

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"), target_waves_Hs; force=true)

        src_bathy = nothing
        target_waves_Hs = nothing
        GC.gc()
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"))
        # Need MPA bathy data as a template for wave data due to unknown crs/format
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

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"), target_waves_Tp; force=true)

        src_bathy = nothing
        target_waves_Tp = nothing
        GC.gc()
    end
end