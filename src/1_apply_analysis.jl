using Rasters, YAXArrays
using GLMakie, GeoMakie

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG

using ImageFiltering
using ImageContrastAdjustment: adjust_histogram, LinearStretching
using ImageMorphology: label_components, component_centroids

using Statistics, StatsBase
using Glob
using ProgressMeter


"""
    suitability_func(threshold::Float64)::Function

Generate function that identifies whether a pixel that has an area that meets the
suitability threshold.
"""
function suitability_func(threshold::Float64)::Function
    function is_suitable(subsection::AbstractMatrix)::Int32
        total = sum(subsection)
        if total == 0.0
            return 0.0
        end

        return Int32((total / length(subsection)) .> threshold)
    end

    return is_suitable
end

# Criteria/suitability functions
function depth_criteria(depth)
    return (-3 .< depth .< 0.0)
end
function slope_criteria(slope)
    return (0.0 .< slope .< 40.0)
end
function supports_coral(p)
    # Manually extracted from Raster Attribute Table
    # 0x00 = 0 = missing value
    # 0x0b = 11 = Sand
    # 0x0c = 12 = Rubble
    # 0x0d = 13 = Rock
    # 0x0f = 15 = Coral/Algae
    return (p .∈ [[0x0d, 0x0f]])
end

# Image filtering functions
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


# Script start
DATA_DIR = "C:/Users/tiwanaga/development/ADRIA_data/spatial_datasets/Bathy data 10m"
RESULT_DIR = "../figs/"

gbr_benthic_path = "$(DATA_DIR)/Benthic/GBR10 GBRMP Benthic.tif"
gbr_benthic = Raster(gbr_benthic_path, crs=EPSG(4326), mappedcrs=EPSG(32755), lazy=true)

regions = [
    "Townsville-Whitsunday",
    "Cairns-Cooktown",
    "Mackay-Capricorn",
    "FarNorthern",
]


@showprogress desc="Analyzing..." for reg in regions
    # Load bathymetry raster
    src_bathy_path = first(glob("*.tif", joinpath(DATA_DIR, "bathy", reg)))
    src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    # Load slope raster
    src_slope_path = first(glob("*.tif", joinpath(DATA_DIR, "slope", reg)))
    src_slope = Raster(src_slope_path, mappedcrs=EPSG(4326), lazy=true)

    # Convert EPSG:32755 to EPSG:4326 simply to crop the larger GBR-wide dataset down
    # to the target region (which is in EPSG:4326).
    # This is to cut down the amount of data to process
    dummy = resample(src_slope; crs=EPSG(4326))
    src_benthic = crop(gbr_benthic; to=dummy)
    dummy = nothing  # mark for garbage collector
    src_benthic = resample(src_benthic, to=src_slope)

    # Source image is of 10m^2 pixels
    # A hectare is 100x100 meters, so we're looking for contiguous areas where at least 95%
    # some proportion (here 75 or 95%) meet criteria of
    # (depth < 3, slope < 40, and habitat is Rock or Coral/Algae).
    suitable = read(
        depth_criteria(src_bathy) .&
        slope_criteria(src_slope) .&
        supports_coral(src_benthic)
    )

    res75 = mapwindow(meets_suitability_threshold(0.75), suitable, (-4:5, -4:5), border=Fill(0)) .|> Gray
    result_raster = res75
    write(joinpath(RESULT_DIR, "$(reg)_suitability75.tiff"), result_raster, force=true)
    res75 = nothing  # mark for garbage collector

    res95 = mapwindow(meets_suitability_threshold(0.95), suitable, (-4:5, -4:5), border=Fill(0)) .|> Gray
    write(joinpath(RESULT_DIR, "$(reg)_suitability95.tiff"), res95, force=true)

    result_raster .= label_components(res95)
    write(joinpath(RESULT_DIR, "$(reg)_grouped_suitability.tiff"), result_raster, force=true)

    # Apply mode filter on ~double the area
    mode_res = mapwindow(c_mode, result_raster, (-9:10, -9:10))
    result_raster .= label_components(mode_res)
    write(joinpath(RESULT_DIR, "$(reg)_grouped_filtered_suitability.tiff"), result_raster, force=true)
end
