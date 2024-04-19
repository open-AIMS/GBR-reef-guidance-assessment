"""Prepare data for analysis."""

using Rasters
using NCDatasets
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

@showprogress dt=10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
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

    if !isfile(joinpath(RESULT_DIR, "$(reg)_waves.tif"))
        if reg == "Townsville-Whitsunday"
            tw_waves_1 = Rasters.Raster(joinpath(DATA_DIR,"waves/TownsvilleWhitsunday_VarWind_ubed_90.nc"); key=:ubed90, lazy=true, missingval=missingval(src_bathy))
            tw_waves_2 = Rasters.Raster(joinpath(DATA_DIR,"waves/TownsvilleWhitsunday_VarWind_ubed_90_additional.nc")) # contains additional information for reefs in the following cropped region
            tw_waves_2_cropped = tw_waves_2[X(Between(675150,699000)), Y(Between(7765000,7755000))]
            tw_waves = Rasters.mosaic(last, tw_waves_1, tw_waves_2_cropped) # inputting values from the cropped region of the second dataset to the total dataset
            tw_waves = Rasters.resample(tw_waves, to=src_bathy)
            write(joinpath(RESULT_DIR, "$(reg)_waves.tif"), tw_waves; force=true)
        else
            target_waves_path = first(glob("*.tif", joinpath(DATA_DIR, "waves", reg)))
            target_waves = Rasters.Raster(target_waves_path, key=:ubed90, lazy=true, missingval=missingval(src_bathy))
            target_waves = Rasters.resample(target_waves, to=src_bathy)
            write(joinpath(RESULT_DIR, "$(reg)_waves.tif"), target_waves; force=true)
        end
    end
end
