using Distributed
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
using ProgressMeter


include("common.jl")


@everywhere begin
    using Rasters
    using ImageCore: Gray
    using ImageFiltering
    using ImageContrastAdjustment: adjust_histogram, LinearStretching
    using ImageMorphology: label_components, component_centroids

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

    function assess_region(reg)
        # Load bathymetry raster
        src_bathy_path = first(glob("*.tif", joinpath(DATA_DIR, "bathy", reg)))
        src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

        # Load slope raster
        src_slope_path = first(glob("*.tif", joinpath(DATA_DIR, "slope", reg)))
        src_slope = Raster(src_slope_path, mappedcrs=EPSG(4326), lazy=true)

        # Load pre-prepared benthic data
        src_benthic_path = "../figs/$(reg)_benthic.tif"
        src_benthic = Raster(src_benthic_path, lazy=true)

        src_geomorphic_path = "../figs/$(reg)_geomorphic.tif"
        src_geomorphic = Raster(src_geomorphic_path, lazy=true)

        # Source image is of 10m^2 pixels
        # A hectare is 100x100 meters, so we're looking for contiguous areas where
        # some proportional area (here 75% or 95%) meet criteria of
        # (-9 <= depth <= -3, slope < 40, and habitat is Rock or Coral/Algae).
        # 75% is assessed for comparison purposes.
        # suitable = read(
        #     depth_criteria(src_bathy) .& slope_criteria(src_slope) .& supports_coral(src_benthic)
        # )

        # See comment above re suitability functions - use of functions breaks `read()`
        suitable_flats = read(
            (src_geomorphic .∈ [FLAT_IDS]) .&
            (src_benthic .∈ [BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0)
        )

        suitable_slopes = read(
            (src_geomorphic .∈ [SLOPE_IDS]) .&
            (src_benthic .∈ [BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0)
        )

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int32, copy(suitable_flats))
        rebuild(result_raster; missingval=0)

        # res75 = mapwindow(suitability_func(0.75), suitable, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # result_raster .= res75
        # write(joinpath(RESULT_DIR, "$(reg)_suitability75.tif"), result_raster; force=true)

        # # Explicitly call GC to free up memory as we go!
        # res75 = nothing
        # GC.gc()

        # res85 = mapwindow(suitability_func(0.85), suitable, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # result_raster .= res85
        # write(joinpath(RESULT_DIR, "$(reg)_suitability85.tif"), result_raster; force=true)
        # res85 = nothing
        # GC.gc()

        # Assess flats
        res95 = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        result_raster .= res95
        write(joinpath(RESULT_DIR, "$(reg)_suitable_flats_95.tif"), result_raster; force=true)

        result_raster .= label_components(res95)
        write(joinpath(RESULT_DIR, "$(reg)_grouped_flats_95.tif"), result_raster; force=true)

        # Assess slopes
        res95 = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        result_raster .= res95
        write(joinpath(RESULT_DIR, "$(reg)_suitable_slopes_95.tif"), result_raster; force=true)

        result_raster .= label_components(res95)
        write(joinpath(RESULT_DIR, "$(reg)_grouped_slopes_95.tif"), result_raster; force=true)
        res95 = nothing
        GC.gc()

        # Apply mode filter (doesn't work as intended)
        # mode_res = mapwindow(c_mode, res95, (5, 5))
        # result_raster .= label_components(mode_res)
        # write(joinpath(RESULT_DIR, "$(reg)_grouped_filtered_suitability.tif"), result_raster; force=true)
        # mode_res = nothing

        # res75 = nothing
        # res95 = nothing
        # result_raster = nothing
    end
end


@showprogress dt=10 desc="Analyzing..." pmap(assess_region, REGIONS)
