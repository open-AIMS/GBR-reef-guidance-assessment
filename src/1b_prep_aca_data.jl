"""
Prepare data for analysis.

Reproject data from WGS84 to UTM Zone 54 - 56.
Ensure all rasters are the same size/shape for each region of interest.
"""

import GeoDataFrames as GDF
using Rasters       # need without common.jl
using GeoInterface

include("common.jl")

aca_output_dir = joinpath(OUTPUT_DIR, "ACA")

# Load bathymetric and turbidity data for separation into regions

aca_bathy_path = "$(ALLEN_ATLAS_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
aca_bathy = Raster(aca_bathy_path, mappedcrs=EPSG(4326), lazy=true)

aca_turbid_path = "$(ALLEN_ATLAS_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
aca_turbid = Raster(aca_turbid_path, lazy=true)

region_features = GDF.read(region_path)

# Reproject region features to ensure they are consistent with aca_bathy
proj_str = ProjString(AG.toPROJ4(AG.importWKT(crs(aca_bathy).val; order=:compliant)))
region_features.geometry = AG.reproject(region_features.geometry, EPSG(4326), proj_str; order=:trad)
region_features[!, :geometry] = Vector{AG.IGeometry}(AG.forceto.(region_features.geometry, AG.wkbMultiPolygon))

@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    # src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    # aca_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    if !isfile(joinpath(aca_output_dir, "$(reg)_bathy.tif"))
        # Cropping the raster speeds up the trim(mask()) step and reduces memory use - checked output in qgis
        target_bathy = Rasters.crop(aca_bathy; to=region_features[reg_idx, :])
        target_bathy = Rasters.trim(mask(target_bathy; with=region_features[reg_idx, :]))

        write(joinpath(aca_output_dir, "$(reg)_bathy.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(aca_output_dir, "$(reg)_turbid.tif"))
        # Using aca_bathy as the basis for resampling for ACA data
        src_aca_bathy = Raster(joinpath(aca_output_dir, "$(reg)_bathy.tif"), lazy=true)

        target_turbid = Rasters.crop(aca_turbid; to = region_features[reg_idx, :])
        target_turbid = Rasters.trim(mask(target_turbid; with=region_features[reg_idx, :]))
        target_turbid = resample(target_turbid, to=src_aca_bathy)

        write(joinpath(aca_output_dir, "$(reg)_turbid.tif"), target_turbid; force=true)

        target_turbid = nothing
        src_aca_bathy = nothing
        GC.gc()
    end
end
