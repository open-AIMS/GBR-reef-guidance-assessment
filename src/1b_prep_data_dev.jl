"""Prepare data for analysis."""

using Rasters
# using NCDatasets
import GeoDataFrames as GDF
import ArchGDAL as AG

using Statistics, StatsBase
using Glob

include("common.jl")

const ALLEN_DATA_DIR = joinpath(DATA_DIR, "..", "AllenAtlas_GBR-20231118074407")

# Reproject benthic data into common CRSs
# Avoids the need to repeatedly reproject within each step

allen_bathy_path = "$(ALLEN_DATA_DIR)/Bathymetry---composite-depth/bathymetry_0.tif"
allen_bathy = Raster(allen_bathy_path, crs=EPSG(4326), lazy=true)

allen_turbid_path = "$(ALLEN_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
allen_turbid = Raster(allen_turbid_path, crs=EPSG(4326), lazy=true)

# using GeoArrays
# using GeoInterface
# function geoJSON_to_tiff(geojson_file::String, tiff_file::String)
#     # Read GeoJSON into a GeoDataFrame
#     gdf = GDF.read(geojson_file)
#     ####
#     # DOESN't WORK BELOW HERE
#     ####
#     # Define the extent and resolution of the raster
#     extent = extrema(gdf) #GeoInterface.get_bounds(gdf)
#     resolution = 0.01  # Set the desired resolution

#     # Convert GeoDataFrame to raster
#     raster = rasterize(gdf, extent, resolution)

#     # Write raster to a TIFF file
#     write(tiff_file, raster)
# end

# allen_benthic_path = "$(ALLEN_DATA_DIR)/Benthic-Map/benthic.geojson"
# # TO DO: convert geojson to raster
# # TO DO: save raster to file
# allen_benthic = Raster(allen_benthic_path, crs=EPSG(4326), lazy=true)

# allen_morphic_path = "$(ALLEN_DATA_DIR)/Geomorphic-Map/geomorphic.geojson"
# # TO DO: convert geojson to raster
# # TO DO: save raster to file
# allen_geomorphic = Raster(allen_morphic_path, crs=EPSG(4326), lazy=true)


# Get polygon of management areas
region_path = joinpath(
    DATA_DIR,
    "zones",
    "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
)
region_features = GDF.read(region_path)

@showprogress dt=10 "Prepping benthic/geomorphic data..." for reg in REGIONS
    reg_idx = occursin.(reg[1:3], region_features.AREA_DESCR)                   # Fn to mask/ignore FarNorthern

    if !isfile(joinpath(RESULT_DIR, "$(reg)_benthic_allen.tif"))
        # Only recreate files if needed
        target_bathy = Rasters.trim(mask(allen_bathy; with=region_features[reg_idx, :]))
        write(joinpath(RESULT_DIR, "$(reg)_benthic_allen.tif"), target_bathy; force=true)

        target_bathy = nothing
        GC.gc()
    end

    if !isfile(joinpath(RESULT_DIR, "$(reg)_turbidity_allen.tif"))
        target_turbidity = Rasters.trim(mask(allen_turbidity; with=region_features[reg_idx, :]))
        write(joinpath(RESULT_DIR, "$(reg)_turbidity_allen.tif"), target_turbidity; force=true)

        target_turbidity = nothing
        GC.gc()
    end

    # if !isfile(joinpath(RESULT_DIR, "$(reg)_benthic.tif"))
    #     # Only recreate files if needed
    #     target_benthic = Rasters.trim(mask(gbr_benthic; with=region_features[reg_idx, :]))
    #     target_benthic = resample(target_benthic, to=src_bathy)
    #     write(joinpath(RESULT_DIR, "$(reg)_benthic.tif"), target_benthic; force=true)

    #     target_benthic = nothing
    #     GC.gc()
    # end
    # if !isfile(joinpath(RESULT_DIR, "$(reg)_geomorphic.tif"))
    #     target_geomorphic = Rasters.trim(mask(gbr_geomorphic; with=region_features[reg_idx, :]))
    #     target_geomorphic = resample(target_geomorphic, to=src_bathy)
    #     write(joinpath(RESULT_DIR, "$(reg)_geomorphic.tif"), target_geomorphic; force=true)

    #     target_geomorphic= nothing
    #     GC.gc()
    # end
end
