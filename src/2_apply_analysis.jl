"""Identify suitable locations for each region."""


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

    function assess_region(reg)
        # Load bathymetry raster
        src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

        # Load slope raster
        src_slope_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        src_slope = Raster(src_slope_path, mappedcrs=EPSG(4326), lazy=true)

        # Load pre-prepared benthic data
        src_benthic_path = joinpath(OUTPUT_DIR, "$(reg)_benthic.tif")
        src_benthic = Raster(src_benthic_path, lazy=true)

        src_geomorphic_path = joinpath(OUTPUT_DIR, "$(reg)_geomorphic.tif")
        src_geomorphic = Raster(src_geomorphic_path, lazy=true)

        src_waves_path = joinpath(OUTPUT_DIR, "$(reg)_waves.tif")
        src_waves = Raster(src_waves_path, lazy=true, crs=crs(src_bathy))

        # Source image is of 10m^2 pixels
        # A hectare is 100x100 meters, so we're looking for contiguous areas where
        # some proportional area (here 75% or 95%) meet criteria of
        # (-9 <= depth <= -3, slope < 40, and habitat is Rock or Coral/Algae).
        # 75% is assessed for comparison purposes.
        # suitable = read(
        #     depth_criteria(src_bathy) .& slope_criteria(src_slope) .& supports_coral(src_benthic)
        # )

        # See comment above re suitability functions - use of functions breaks `read()`
        # Assess flats
        suitable_flats = read(
            (src_geomorphic .∈ [FLAT_IDS]) .&
            (src_benthic .∈ [BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0) .&
            (0.0 .<= src_waves .<= 1.0)
        )

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_flats))
        rebuild(result_raster; missingval=0)

        # # 85% threshold
        # res = mapwindow(suitability_func(0.85), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_flats_85.tif")
        # _write_data(fpath, res, result_raster)

        # fpath = joinpath(OUTPUT_DIR, "$(reg)_grouped_flats_85.tif")
        # _write_data(fpath, res, result_raster)

        # 95% threshold
        res = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_flats_95.tif")
        _write_data(fpath, res, result_raster)

        fpath = joinpath(OUTPUT_DIR, "$(reg)_grouped_flats_95.tif")
        _write_data(fpath, res, result_raster)

        # Assess slopes
        suitable_slopes = read(
            (src_geomorphic .∈ [SLOPE_IDS]) .&
            (src_benthic .∈ [BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0) .&
            (0.0 .<= src_waves .<= 1.0)
        )

        # # 85% threshold
        # res = mapwindow(suitability_func(0.85), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_slopes_85.tif")
        # _write_data(fpath, res, result_raster)

        # fpath = joinpath(OUTPUT_DIR, "$(reg)_grouped_slopes_85.tif")
        # _write_data(fpath, res, result_raster)

        # 95% threshold
        res = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_slopes_95.tif")
        _write_data(fpath, res, result_raster)

        fpath = joinpath(OUTPUT_DIR, "$(reg)_grouped_slopes_95.tif")
        _write_data(fpath, res, result_raster)

        res = nothing
        GC.gc()
    end
end


@showprogress dt = 10 desc = "Analyzing..." pmap(assess_region, REGIONS)
