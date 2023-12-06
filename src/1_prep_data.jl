using Rasters
import GeoDataFrames as GDF
import ArchGDAL as AG

using Statistics, StatsBase
using Glob

include("common.jl")


# Reproject benthic data into common CRSs
# Avoids the need to repeatedly reproject within each step

gbr_benthic_path = "$(DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
gbr_benthic = Raster(gbr_benthic_path, crs=EPSG(4326), lazy=true)


gbr_morphic_path = "$(DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
gbr_geomorphic = Raster(gbr_morphic_path, crs=EPSG(4326), lazy=true)


# Get polygon of management areas
region_path = joinpath(
    DATA_DIR,
    "zones",
    "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
)
region_features = GDF.read(region_path)

@showprogress dt=10 "Prepping benthic/geomorphic data..." for reg in REGIONS
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)

    src_bathy_path = first(glob("*.tif", joinpath(DATA_DIR, "bathy", reg)))
    src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    if !isfile(joinpath(RESULT_DIR, "$(reg)_benthic.tif"))
        # Only recreate files if needed
        target_benthic = Rasters.trim(mask(gbr_benthic; with=region_features[reg_idx, :]))
        target_benthic = resample(target_benthic, to=src_bathy)
        write(joinpath(RESULT_DIR, "$(reg)_benthic.tif"), target_benthic; force=true)

        target_benthic = nothing
        GC.gc()
    end

    if !isfile(joinpath(RESULT_DIR, "$(reg)_geomorphic.tif"))
        target_geomorphic = Rasters.trim(mask(gbr_geomorphic; with=region_features[reg_idx, :]))
        target_geomorphic = resample(target_geomorphic, to=src_bathy)
        write(joinpath(RESULT_DIR, "$(reg)_geomorphic.tif"), target_geomorphic; force=true)

        target_geomorphic= nothing
        GC.gc()
    end
end
