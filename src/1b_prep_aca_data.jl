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

# Reproject benthic data into common CRSs
# Avoids the need to repeatedly reproject within each step

aca_bathy_path = "$(ALLEN_ATLAS_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
aca_bathy = Raster(aca_bathy_path, mappedcrs=EPSG(4326), lazy=true)

aca_turbid_path = "$(ALLEN_ATLAS_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
aca_turbid = Raster(aca_turbid_path, crs=EPSG(4326), lazy=true)

region_features = GDF.read(region_path)

@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    # src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    # aca_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    if !isfile(joinpath(aca_output_dir, "$(reg)_bathy.tif"))
        target_bathy = Rasters.trim(mask(aca_bathy; with=region_features[reg_idx, :]))
        # target_bathy = resample(target_bathy, to=aca_bathy)
        write(joinpath(aca_output_dir, "$(reg)_bathy.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(OUTPUT_DIR, "$(reg)_turbid.tif"))
        target_turbid = Rasters.trim(mask(aca_turbid; with=region_features[reg_idx, :]))
        # target_turbid = resample(target_turbid, to=aca_bathy)
        write(joinpath(OUTPUT_DIR, "$(reg)_turbid.tif"), target_turbid; force=true)

        target_turbid = nothing
        GC.gc()
    end
end
