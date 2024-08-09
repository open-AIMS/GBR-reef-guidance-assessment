"""
Prepare data for analysis by processing ACA files for each GBRMPA management region.

Filter desired FLAT, SLOPE and BENTHIC IDs to save time in 2b analysis.
Crop GBR-wide ACA rasters into management regions.
Reproject all data from WGS84 / UTM Zone 54 - 56 into consistent CRS (GDA2020).
Ensure all rasters are the same size/shape for each region of interest with the same
values used to indicate no data.

The general approach is to crop down to an extent of the region, then trim/mask to just the
areas of interest. Processing data in this way minimized the amount of data being handled.

Still, the largest chunk of time is spent resampling to the same projections and writing
data to file (as it takes time to compress the data). The datasets can be in the GBs without
compression.
"""

include("common.jl")

# Loading regions_4326 for cropping of vector and raster data.
regions_4326 = GDF.read(REGION_PATH_4326)

# 1. Processing of geojson files into smaller-GBRMPA regions

# Loading GBR-wide polygon data
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

    flats_fn = joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg).gpkg")
    if !isfile(flats_fn)
        flat_is_in_region = AG.contains.([region_4326_geom], target_flat_poly.geometry)
        target_flats_reg = target_flat_poly[flat_is_in_region, :]

        target_flats_reg.geometry = AG.reproject(
            target_flats_reg.geometry,
            crs(region_4326_geom),
            GDA2020_crs;
            order=:trad
        )

        GDF.write(flats_fn, target_flats_reg; crs=EPSG(7844))
    end

    slopes_fn = joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg).gpkg")
    if !isfile(slopes_fn)
        slope_is_in_region = AG.contains.([region_4326_geom], target_slope_poly.geometry)
        target_slopes_reg = target_slope_poly[slope_is_in_region, :]

        target_slopes_reg.geometry = AG.reproject(
            target_slopes_reg.geometry,
            crs(region_4326_geom),
            GDA2020_crs;
            order=:trad
        )

        GDF.write(slopes_fn, target_slopes_reg; crs=EPSG(7844))
    end

    benthic_fn = joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg).gpkg")
    if !isfile(benthic_fn)
        ben_is_in_region = AG.contains.([region_4326_geom], benthic_poly.geometry)
        benthic_reg = benthic_poly[ben_is_in_region, :]

        benthic_reg.geometry = AG.reproject(
            benthic_reg.geometry,
            crs(region_4326_geom),
            GDA2020_crs;
            order=:trad
        )

        GDF.write(benthic_fn, benthic_reg; crs=EPSG(7844))
    end
end

# 2. Processing of reef-wide rasters into management region scale rasters, and reprojection
# of bathymetry, turbidity and wave raster data.

# Loading GBR-wide data
aca_bathy_path = "$(ACA_DATA_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
aca_bathy = Raster(aca_bathy_path, mappedcrs=EPSG(4326), lazy=true)

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping bathymetry/turbidity/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    base_bathy_fn = joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif")
    if !isfile(base_bathy_fn)
        target_bathy = Rasters.crop(aca_bathy; to=regions_4326[reg_idx_4326, :])
        target_bathy = Rasters.trim(mask(target_bathy; with=regions_4326[reg_idx_4326, :]))
        target_bathy = set_consistent_missing(target_bathy, -9999.0)

        write(
            base_bathy_fn,
            resample(
                target_bathy;
                crs=GDA2020_crs
            )
        )

        target_bathy = nothing
        force_gc_cleanup()
    end
end
