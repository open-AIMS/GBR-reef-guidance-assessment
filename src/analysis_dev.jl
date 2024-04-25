using Rasters, YAXArrays

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG

using GLMakie, GeoMakie

using ImageFiltering
using ImageContrastAdjustment: adjust_histogram, LinearStretching
using ImageMorphology: label_components, component_centroids
using Statistics, StatsBase


# Examine Bathymetry data (10m resolution)
MPA_DATA_DIR = "C:/Users/tiwanaga/development/ADRIA_data/spatial_datasets/Bathy data 10m"

# tsv_data = joinpath(MPA_DATA_DIR, "Townsville-Whitsunday/SDB_AUS_TownsvilleWhitsunday_EOMAP_20180603_20180903_10m_MSL_geotiff.tif")
cairns_data = joinpath(MPA_DATA_DIR, "bathy/Cairns-Cooktown/SDB_AUS_CairnsCooktown_EOMAP_20180514_20181016_10m_MSL_geotiff.tif")
cairns_slope_data = joinpath(MPA_DATA_DIR, "slope/Cairns-Cooktown/slope_cnscook.tif")

moore_dpkg_path = "C:/Users/tiwanaga/development/ADRIA_data/Moore_2023-09-14_rc/site_data/Moore_2023-09-14_rc.gpkg"
gbr_dpkg_path = "C:/Users/tiwanaga/development/RME/rme_ml_2023_03_30b/data_files/region/reefmod_gbr.gpkg"

cairns_bathy = Raster(cairns_data, mappedcrs=EPSG(4326), lazy=true)
cairns_slope = Raster(cairns_slope_data, mappedcrs=EPSG(4326), lazy=true)
# moore_gdf = GDF.read(moore_dpkg)

gbr_benthic_path = "C:/Users/tiwanaga/development/ADRIA_data/spatial_datasets/Bathy data 10m/benthic/GBR10 GBRMP Benthic.tif"
gbr_benthic = Raster(gbr_benthic_path, crs=EPSG(4326), mappedcrs=EPSG(32755), lazy=true)

# Convert EPSG:32755 to EPSG:4326 simply to crop the larger GBR-wide dataset down
# This is to cut down the amount of data to process
dummy = resample(cairns_slope; crs=EPSG(4326))
cairns_benthic = crop(gbr_benthic; to=dummy)
dummy = nothing
cairns_benthic = resample(cairns_benthic, to=cairns_slope)


"""Identify pixel that has an area that meets the suitability threshold"""
function meets_suitability_threshold(threshold::Float64)::Function
    function is_suitable(subsection::AbstractMatrix)::Int32
        total = sum(subsection)
        if total == 0.0
            return 0.0
        end

        return Int32((total / length(subsection)) .> threshold)
    end

    return is_suitable
end


gbr_dpkg = GDF.read(gbr_dpkg_path)
gbr_dpkg.geom = GDF.reproject(gbr_dpkg.geom, GFT.EPSG(4326), GFT.EPSG(32755), order=:trad)

# replace_missing!(cairns_slope, 9999)

# Manually extracted from Raster Attribute Table
# 0x00 = 0 = missing value
# 0x0b = 11 = Sand
# 0x0c = 12 = Rubble
# 0x0d = 13 = Rock
# 0x0f = 15 = Coral/Algae

# Source image is of 10m^2 pixels
# A hectare is 100x100, so we're looking for a contiguous area of at least 10x10 pixels
# that fit our criteria (depth < 3, slope < 40, Rock or Coral/Algae).
cairns_1_suitability = read((-3 .< cairns_bathy .< 0.0) .& (0.0 .< cairns_slope .< 40.0) .& (cairns_benthic .∈ [[0x0d, 0x0f]]))

# If specifying window size with integers, the template/window must be odd as it will be
# centered.
# (3, 3) becomes -1:1 (i.e., -1, 0, 1)
# We can also specify an interval like -4:5
# which would something like:
# 1, 1, 1, 1, 0, 1, 1, 1, 1, 1
# where 0 is the center of the window (only used here as an example,
# it is not given a 0 value in the analysis).
@info "Window filtering"
res75 = mapwindow(meets_suitability_threshold(0.75), cairns_1_suitability, (-4:5, -4:5), border=Fill(0)) .|> Gray
result_raster = convert.(Int32, copy(cairns_1_suitability))
result_raster .= res75
write("cairns_suitability75.tiff", result_raster, force=true)

res75 = nothing

res95 = mapwindow(meets_suitability_threshold(0.95), cairns_1_suitability, (-4:5, -4:5), border=Fill(0)) .|> Gray
write("cairns_suitability95.tiff", res95, force=true)

result_raster .= label_components(res95)
write("cairns_grouped_suitability.tiff", result_raster, force=true)

centroids = component_centroids(result_raster)[2:end]

function c_median(x)
    if sum(x) == 0
        return 0
    end

    return floor(Int64, median(x))
end

function c_mode(x)
    if sum(x) == 0
        return 0
    end

    return mode(x)
end

# Apply mode filter
imgm = mapwindow(c_mode, result_raster, (-9:10, -9:10))

result_raster .= label_components(imgm)

write("cairns_grouped_filtered_suitability.tiff", result_raster, force=true)
