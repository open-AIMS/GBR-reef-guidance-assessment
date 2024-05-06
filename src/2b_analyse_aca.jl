"""Identify suitable locations for each region."""

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

    function analyze_allen()
        bathy_rst = Raster(joinpath(ACA_DATA_DIR, "Bathymetry---composite-depth", "bathymetry_0.tif"), lazy=true)
        turbid_rst = Raster(joinpath(ACA_DATA_DIR, "Turbidity-Q3-2023", "turbidity-quarterly_0.tif"), lazy=true) # mappedcrs=EPSG(4326), 

        suitable_raster = read(      # this uses a lot of memory 
            (200.0 .<= bathy_rst .<= 900.0)  # depth criteria
            .& (turbid_rst .<= 52)    # LOW Turbidity = “52” in the turbidity maps = “5.2 FNU”
        )    
        convert.(Int16, suitable_raster)
        rebuild(suitable_raster; missingval=0)      
        
        # # Source image is of 10m^2 pixels
        # # A hectare is 100x100 meters, so we're looking for contiguous areas where
        # # some proportional area (here 75% or 95%) meet criteria of
        # # (-9 <= depth <= -3, slope < 40, and habitat is Rock or Coral/Algae).
        # # 75% is assessed for comparison purposes.
        # # suitable = read(
        # #     depth_criteria(src_bathy) .& slope_criteria(src_slope) .& supports_coral(src_benthic)
        # # )

        # # See comment above re suitability functions - use of functions breaks `read()`

        # mask suitable areas by raster files with geomorphic and benthic data previously saved as gpkg
        benthic_morphic_poly = GDF.read(first(glob("*.gpkg", ACA_OUTPUT_DIR)))
        suitable_areas = Rasters.mask(suitable_raster; benthic_morphic_poly) .|> Gray

        #####################
        geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson")) # first(glob("*.geojson", joinpath(ACA_DATA_DIR, "Geomorphic-Map")))
        # flats
        # target geomorphic zones by flat and slope
        target_flats = geomorphic_poly[geomorphic_poly.class .∈ [ACA_FLAT_IDS], :]      # Deep lagoon, Shallow Lagoon, Terrestrial Reef flat, Inner Reef Flat, Outer Reef Flat, Plateau
        suitable_flats = Rasters.mask(suitable_areas; with=target_flats) .|> Gray
        # Need a copy of raster data type to support writing to `tif`
        result_raster_flats = convert.(Int16, copy(suitable_flats))
        rebuild(result_raster_flats; missingval=0)


        # # 85% threshold
        # res = mapwindow(suitability_func(0.85), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_flats_85.tif")
        # _write_data(fpath, res, result_raster)

        # fpath = joinpath(OUTPUT_DIR, "$(reg)_grouped_flats_85.tif")
        # _write_data(fpath, res, result_raster)

        # 95% threshold
        res = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_flats_95.tif")
        _write_data(fpath, res, result_raster)

        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_grouped_flats_95.tif")
        _write_data(fpath, res, result_raster)

        suitable_flats = nothing
        GC.gc()     # delete?

        res = nothing
        GC.gc()

        #####################
        # slopes
        target_slopes = geomorphic_poly[geomorphic_poly.class .∈ [ACA_SLOPE_IDS], :]    # Sheltered Reef Slope, Reef Slope, Back Reef Slope
        suitable_slopes = Rasters.mask(suitable_areas; with=target_slopes) .|> Gray
        # Need a copy of raster data type to support writing to `tif`
        result_raster_slopes = convert.(Int16, copy(suitable_slopes))
        rebuild(result_raster_slopes; missingval=0)

        geomorphic_poly = nothing
        GC.gc()

        # # 85% threshold
        # res = mapwindow(suitability_func(0.85), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_slopes_85.tif")
        # _write_data(fpath, res, result_raster)

        # fpath = joinpath(OUTPUT_DIR, "$(reg)_grouped_slopes_85.tif")
        # _write_data(fpath, res, result_raster)

        # 95% threshold
        res = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes_95.tif")
        _write_data(fpath, res, result_raster)

        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_grouped_slopes_95.tif")
        _write_data(fpath, res, result_raster)

        suitable_slopes = nothing
        GC.gc()     # delete?

        res = nothing
        GC.gc()
    end
end
    
@showprogress dt = 10 desc = "Analyzing..." pmap(analyze_allen)