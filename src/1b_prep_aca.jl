"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

include("common.jl")

# 1. Processing of geojson files into a reef-wide geopackage
if !isfile(first(glob("*.gpkg", ACA_OUTPUT_DIR)))
    geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson"))
    benthic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Benthic-Map", "benthic.geojson"))
    reef_poly = GDF.read(joinpath(ACA_DATA_DIR, "Reef-Extent", "reefextent.geojson"))

    target_geomorphic_features = geomorphic_poly.class .∈ Ref([ACA_FLAT_IDS..., ACA_SLOPE_IDS])
    target_benthic_features = benthic_poly.class .∈ Ref(ACA_BENTHIC_IDS)

    geomorphic_poly = geomorphic_poly[target_geomorphic_features, :]
    benthic_poly = benthic_poly[target_benthic_features, :]
    GC.gc()

    # Assuming we skip any reefs where we do not have data across all criteria
    skipped_reefs = fill("", size(reef_poly, 1))
    @floop for (idx, reef) in enumerate(eachrow(reef_poly))
        reef_geom = reef.geometry
        has_geomorphic = AG.intersects.([reef_geom], geomorphic_poly.geometry)
        if all(has_geomorphic .== false)
            skipped_reefs[idx] = "No geomorphic"
            continue
        end

        has_benthic = AG.intersects.([reef_geom], benthic_poly.geometry)
        if all(has_benthic .== false)
            skipped_reefs[idx] = "No benthic"
            continue
        end
    end
    GDF.write(joinpath(ACA_OUTPUT_DIR, "aca_benthic_geomorphic.gpkg"), reef_poly[skipped_reefs.=="", :])
end

# 2. Processing of reef-wide rasters into smaller-GBRMPA regions
aca_bathy_path = "$(ACA_DATA_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
aca_bathy = Raster(aca_bathy_path, mappedcrs=EPSG(4326), lazy=true)

aca_turbid_path = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
aca_turbid = Raster(aca_turbid_path, lazy=true)

region_features = GDF.read(REGION_PATH)

# Reproject region features to ensure they are consistent with aca_bathy
proj_str = ProjString(AG.toPROJ4(AG.importWKT(crs(aca_bathy).val; order=:compliant)))
region_features.geometry = AG.reproject(region_features.geometry, EPSG(4326), proj_str; order=:trad)
region_features[!, :geometry] = Vector{AG.IGeometry}(AG.forceto.(region_features.geometry, AG.wkbMultiPolygon))

@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
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
end
