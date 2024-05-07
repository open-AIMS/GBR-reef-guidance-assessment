"""
Identify suitable locations based on bathymetric, turbidity, geomorphic and benthic criteria.
Use raster and polygon data to analyze suitabile areas.
Output raster files with areas meeting specified suitability thresholds.
"""

using Rasters

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG

using ImageCore: Gray
using ImageFiltering
using ImageContrastAdjustment: adjust_histogram, LinearStretching
using ImageMorphology: label_components, component_centroids

using Statistics, StatsBase
using Glob


include("common.jl")

@everywhere begin

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

            return Int32((total / length(subsection)) .>= threshold)
        end

        return is_suitable
    end
    function prop_suitable(subsection::AbstractMatrix)::Float32
        total = sum(subsection)
        if total == 0.0
            return 0.0
        end
        return Float32((total / length(subsection)))
    end

    function _write_data(fpath::String, data, cache)::Nothing
        if !isfile(fpath)
            if occursin("grouped_", fpath)
                cache .= label_components(data)
            else
                cache .= data
            end
            write(fpath, cache; force=true)
        end
        return nothing
    end

    # Criteria/suitability functions
    # For some reason, wrapping conditions in functions causes `read()` to fail
    # Kept here in case it becomes useful later.
    function depth_criteria(depth)
        return (-9 .<= depth .<= -2.0)
    end
    function slope_criteria(slope)
        return (0.0 .< slope .< 40.0)
    end
    function supports_coral(benthic)
        # Manually extracted from Raster Attribute Table
        # 0x00 = 0 = missing value
        # 0x0b = 11 = Sand
        # 0x0c = 12 = Rubble
        # 0x0d = 13 = Rock
        # 0x0f = 15 = Coral/Algae
        return (benthic .∈ [[0x0d, 0x0f]])
    end

    # Image filtering functions
    function c_median(x)
        if sum(x) == 0
            return 0
        end

        return floor(Int64, median(x))
    end

    function c_mode(x)
        # Check center of odd-sized square matrix
        d = size(x)
        centroid = ceil.(Int64, d .* 0.5)
        if x[centroid...] == 0
            return 0
        end

        total = sum(x)
        if total <= 3
            return 0
        end

        return mode(x)
    end
end

function analyze_allen()
    bathy_rst = Raster(joinpath(ACA_DATA_DIR, "Bathymetry---composite-depth", "bathymetry_0.tif"), lazy=true)
    turbid_rst = Raster(joinpath(ACA_DATA_DIR, "Turbidity-Q3-2023", "turbidity-quarterly_0.tif"), lazy=true) # mappedcrs=EPSG(4326), 

    # # Source image is of 10m^2 pixels
    # # A hectare is 100x100 meters, so we're looking for contiguous areas where
    # # some proportional area (here 75% or 95%) meet criteria
    # # 200 <= depth <= 900 cm, and turbidty is no more than LOW.

    # # See comment above re suitability functions - use of functions breaks `read()`
    
    suitable_raster = read(
        (200.0 .<= bathy_rst .<= 900.0)
        .& (turbid_rst .<= 52)
    )    
    potential_bathy_turbid = convert.(UInt8, suitable_raster)
    suitable_raster = nothing
    GC.gc()
    rebuild(potential_bathy_turbid; missingval=0)        

    # mask suitable areas by raster files with geomorphic and benthic data previously saved as gpkg
    benthic_morphic_poly = GDF.read(first(glob("*.gpkg", ACA_OUTPUT_DIR)))
    suitable_areas = Rasters.trim(mask(potential_bathy_turbid; with=benthic_morphic_poly))

    benthic_morphic_poly = nothing
    GC.gc()

    # Assess geomorphic suitability
    geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson"))
    
    # flats
    target_flats = geomorphic_poly[geomorphic_poly.class .∈ [ACA_FLAT_IDS], :]
    suitable_flats = Rasters.mask(suitable_areas; with=target_flats)
    # Need a copy of raster data type to support writing to `tif`
    result_raster_flats = convert.(UInt8, copy(suitable_flats))
    rebuild(result_raster_flats; missingval=0)

    # # 85% threshold
    # res = mapwindow(suitability_func(0.85), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
    # fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_flats_85.tif")
    # _write_data(fpath, res, result_raster)

    # 95% threshold
    res = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
    fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_flats_95.tif")
    _write_data(fpath, res, result_raster)

    suitable_flats = nothing
    GC.gc()     # delete?

    res = nothing
    GC.gc()

    #####################
    # slopes
    target_slopes = geomorphic_poly[geomorphic_poly.class .∈ [ACA_SLOPE_IDS], :]    # Sheltered Reef Slope, Reef Slope, Back Reef Slope
    suitable_slopes = Rasters.mask(suitable_areas; with=target_slopes)
    geomorphic_poly = nothing
    suitable_areas = nothing
    GC.gc()

    # Need a copy of raster data type to support writing to `tif`
    result_raster_slopes = convert.(UInt8, copy(suitable_slopes))
    rebuild(result_raster_slopes; missingval=0)

    # # 85% threshold
    # res = mapwindow(suitability_func(0.85), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
    # fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_slopes_85.tif")
    # _write_data(fpath, res, result_raster)

    # 95% threshold
    res = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
    fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_slopes_95.tif")
    _write_data(fpath, res, result_raster)

    suitable_slopes = nothing
    GC.gc()     # delete?

    res = nothing
    GC.gc()
end
    
# @showprogress dt = 10 desc = "Analyzing..." pmap(analyze_allen)

analyze_allen()