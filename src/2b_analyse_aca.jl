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

aca_output_dir = joinpath(OUTPUT_DIR, "ACA")

@everywhere begin
    """
    suitability_func(threshold::Float64)::Function

    Generate function that identifies whether a pixel that has an area that meets the
    suitability threshold.
    """

    aca_output_dir = joinpath(OUTPUT_DIR, "ACA")

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
        aca_bathy_path = joinpath(aca_output_dir, "$(reg)_bathy.tif")
        allen_bathy = Raster(aca_bathy_path, crs=EPSG(4326), lazy=true)    
    
        # allen_turbid_path = joinpath(aca_output_dir, "$(reg)_turbid.tif")
        # allen_turbi = Raster(allen_turbid_path, crs=EPSG(4326), lazy=true)    
    
        suitable_areas = read(
            (-9.0 .<= allen_bathy .<= -2.0)  # depth criteria
            # .& (allen_turbi .<= 52)    # LOW Turbidity = “52” in the turbidity maps = “5.2 FNU”
        )    
        allen_bathy = nothing
        # allen_turbi = nothing
        GC.gc()    
    
        flat_aca_ids = ["Inner Reef Flat", "Outer Reef Flat", "Plateau"]           # identify flats
        slope_aca_ids = ["Sheltered Reef Slope", "Reef Slope", "Back Reef Slope"]  # identify slopes        allen_bathy_path = first(glob("*.tif", joinpath(ALLEN_ATLAS_DIR, "Bathymetry---composite-depth")))
        
        allen_geo_path = first(glob("*.geojson", joinpath(ALLEN_ATLAS_DIR, "Geomorphic-Map")))
        allen_geo = GDF.read(allen_geo_path)

        all_flats = filter(row -> row[:class] in flat_aca_ids, allen_geo)
        all_slopes = filter(row -> row[:class] in slope_aca_ids, allen_geo)
        allen_geo = nothing
        GC.gc()

        # Attempted fix for Rasters.mask ERROR: InexactError: Bool(-9999.0)
        # suitable_areas = ismissing.(suitable_areas) .| suitable_areas 

        suitable_flats = Rasters.mask(suitable_areas; with=all_flats)
        suitable_slopes = Rasters.mask(suitable_areas; with=all_slopes)

        benthic_aca_ids = ["Rock", "Coral/Algae"]
        aca_benthic_path = first(glob("*.geojson", joinpath(ALLEN_ATLAS_DIR, "Benthic-Map")))
        aca_benthic = GDF.read(aca_benthic_path)

        suitable_benthic = filter(row -> row[:class] in benthic_aca_ids, aca_benthic)

        ###############################################
        #### RUNS TO HERE ####
        ###############################################
    
    
        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_flats))
        rebuild(result_raster; missingval=0)    
    
        ####

        # # Source image is of 10m^2 pixels
        # # A hectare is 100x100 meters, so we're looking for contiguous areas where
        # # some proportional area (here 75% or 95%) meet criteria of
        # # (-9 <= depth <= -3, slope < 40, and habitat is Rock or Coral/Algae).
        # # 75% is assessed for comparison purposes.
        # # suitable = read(
        # #     depth_criteria(src_bathy) .& slope_criteria(src_slope) .& supports_coral(src_benthic)
        # # )

        # # See comment above re suitability functions - use of functions breaks `read()`

        # # Assess flats
        # suitable_flats = read(
        #     (src_geomorphic .∈ [FLAT_IDS]) .&
        #     (src_benthic .∈ [BENTHIC_IDS]) .&
        #     (-9.0 .<= src_bathy .<= -2.0) .&
        #     (0.0 .<= src_slope .<= 40.0) 
        #     #.& (0.0 .<= src_waves .<= 1.0)
        # )

        # # Need a copy of raster data type to support writing to `tif`
        # result_raster = convert.(Int16, copy(suitable_flats))
        # rebuild(result_raster; missingval=0)

        # suitable_flats = nothing
        # GC.gc()

        ####

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

        suitable_slopes = nothing
        GC.gc()

        res = nothing
        GC.gc()
    end
end


@showprogress dt = 10 desc = "Analyzing..." pmap(assess_region, REGIONS)
