"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

include("common.jl")

# 1. Processing of geojson files into smaller-GBRMPA regions
# Loading GBR-wide data
geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson"))
benthic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Benthic-Map", "benthic.geojson"))
reef_poly = GDF.read(joinpath(ACA_DATA_DIR, "Reef-Extent", "reefextent.geojson"))

target_flats = geomorphic_poly.class .∈ Ref(ACA_FLAT_IDS)
target_flat_poly = geomorphic_poly[target_flats, :]

target_slopes = geomorphic_poly.class .∈ Ref(ACA_SLOPE_IDS)
target_slope_poly = geomorphic_poly[target_slopes, :]

target_benthic_features = benthic_poly.class .∈ Ref(ACA_BENTHIC_IDS)
target_benthic_poly = benthic_poly[target_benthic_features, :]

# Subsetting data by region
region_features = GDF.read(REGION_PATH)

@floop for region in eachrow(region_features)
    reg_name = REGIONS[occursin.(region.AREA_DESCR[1:3], REGIONS)][1]
    region_geom = region.geometry

    if !isfile(joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg_name).gpkg"))
        flat_is_in_region = AG.contains.([region], target_flat_poly.geometry)
        target_flats_reg = target_flat_poly[flat_is_in_region, :]

        GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg_name).gpkg"), target_flats_reg)
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg_name).gpkg"))
        slope_is_in_region = AG.contains.([region], target_slope_poly.geometry)
        target_slopes_reg = target_slope_poly[slope_is_in_region, :]

        GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg_name).gpkg"), target_slopes_reg)
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg_name).gpkg"))
        ben_is_in_region = AG.contains.([region], benthic_poly.geometry)
        benthic_reg = benthic_poly[ben_is_in_region, :]

        GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg_name).gpkg"), benthic_reg)
    end
end

# 2. Processing of reef-wide rasters into smaller-GBRMPA regions
# Loading GBR-wide data
aca_bathy_path = "$(ACA_DATA_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
aca_bathy = Raster(aca_bathy_path, mappedcrs=EPSG(4326), lazy=true)

aca_turbid_path = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
aca_turbid = Raster(aca_turbid_path, lazy=true)

# Subsetting data by region
region_features = GDF.read(REGION_PATH)

# Reproject region features to ensure they are consistent with aca_bathy
proj_str = ProjString(AG.toPROJ4(AG.importWKT(crs(aca_bathy).val; order=:compliant)))

region_features.geometry = AG.reproject(region_features.geometry, EPSG(4326), proj_str; order=:trad)
region_features[!, :geometry] = Vector{AG.IGeometry}(AG.forceto.(region_features.geometry, AG.wkbMultiPolygon))

@showprogress dt = 10 "Prepping bathymetric/turbidity/wave data..." for reg in REGIONS
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"))
        # Cropping the raster speeds up the trim(mask()) step and reduces memory use - checked output in qgis
        target_bathy = Rasters.crop(aca_bathy; to=region_features[reg_idx, :])
        target_bathy = Rasters.trim(mask(target_bathy; with=region_features[reg_idx, :]))

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"))
        # Using aca_bathy as the basis for resampling for ACA data
        src_aca_bathy = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"), lazy=true)

        target_turbid = Rasters.crop(aca_turbid; to=region_features[reg_idx, :])
        target_turbid = Rasters.trim(mask(target_turbid; with=region_features[reg_idx, :]))
        target_turbid = resample(target_turbid, to=src_aca_bathy)

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"), target_turbid; force=true)

        target_turbid = nothing
        src_aca_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"))
        target_waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))
        target_waves_Hs = Raster(target_waves_Hs_path, key=:Hs90, crs=crs(src_aca_bathy), mappedcrs=EPSG(4326), lazy=true)

        # Extend bounds of wave data to match bathymetry if needed
        if size(src_aca_bathy) !== size(target_waves_Hs)
            target_waves_Hs = extend(crop(target_waves_Hs; to=src_aca_bathy); to=AG.extent(src_aca_bathy))
            @assert size(src_aca_bathy) == size(target_waves_Hs)

            replace_missing!(target_waves_Hs, -9999.0)
        end

        tmp_Hs = copy(src_aca_bathy)

        # Replace data (important: flip the y-axis!)
        tmp_Hs.data .= coalesce.(target_waves_Hs.data[:, end:-1:1], -9999.0)
        target_waves_Hs = tmp_Hs

        # Set to known missing value
        target_waves_Hs.data[target_waves_Hs.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves_Hs, -9999.0)

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"), target_waves_Hs; force=true)

        target_waves_Hs = nothing
        GC.gc()
    end

    if !isfile(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"))
        target_waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, Tp, reg)))
        target_waves_Tp = Raster(target_waves_Tp_path, key=:Tp90, crs=crs(src_aca_bathy), mappedcrs=EPSG(4326), lazy=true)

        # Extend bounds of wave data to match bathymetry if needed
        if size(src_aca_bathy) !== size(target_waves_Tp)
            target_waves_Tp = extend(crop(target_waves_Tp; to=src_aca_bathy); to=AG.extent(src_aca_bathy))
            @assert size(src_aca_bathy) == size(target_waves_Tp)

            replace_missing!(target_waves_Tp, -9999.0)
        end

        tmp_Tp = copy(src_aca_bathy)

        # Replace data (important: flip the y-axis!)
        tmp_Tp.data .= coalesce.(target_waves_Tp.data[:, end:-1:1], -9999.0)
        target_waves_Tp = tmp_Tp

        # Set to known missing value
        target_waves_Tp.data[target_waves_Tp.data .< -9999.0] .= -9999.0
        replace_missing!(target_waves_Tp, -9999.0)

        write(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"), target_waves_Tp; force=true)

        target_waves_Tp = nothing
        GC.gc()
    end

end
