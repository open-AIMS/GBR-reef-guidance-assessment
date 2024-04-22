"""Identify suitable locations per region."""


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

# using DataFrames
# using Makie

const ALLEN_DATA_DIR = joinpath(DATA_DIR, "..", "AllenAtlas_GBR-20231118074407")

const ALLEN_RESULT_DIR = joinpath(ALLEN_DATA_DIR, "outputs")

const LAYERS =  ["Bathymetry---composite-depth", # "Benthic-Map", "Geomorphic-Map", # these files to be converted to TIFs
                    "Turbidity-Q3-2023"]          #readdir(DATA_DIR)    # Criteria = all files/folders in the directory


                    function mask_by_class(raster::AbstractRaster; geo_df::GeoDataFrames.GeoDataFrame, classes::Vector{String})
    # Filter the GeoDataFrame based on the specified classes
    filtered_geo_df = filter(row -> in(row.class, classes), geo_df)
    
    # Mask the raster with the filtered GeoDataFrame
    masked_raster = Rasters.mask(raster; with=filtered_geo_df)
    
    return masked_raster
end


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

    # # Criteria/suitability functions
    # function assess_region(reg)
    #     # Load bathymetry raster
    #     src_bathy_path = first(glob("*.tif", joinpath(DATA_DIR, "bathy", reg)))
    #     src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    #     # Load slope raster
    #     src_slope_path = first(glob("*.tif", joinpath(DATA_DIR, "slope", reg)))
    #     src_slope = Raster(src_slope_path, mappedcrs=EPSG(4326), lazy=true)

    #     # Load pre-prepared benthic data
    #     src_benthic_path = "../figs/$(reg)_benthic.tif"
    #     src_benthic = Raster(src_benthic_path, lazy=true)

    #     src_geomorphic_path = "../figs/$(reg)_geomorphic.tif"
    #     src_geomorphic = Raster(src_geomorphic_path, lazy=true)

    #     # Source image is of 10m^2 pixels
    #     # A hectare is 100x100 meters, so we're looking for contiguous areas where
    #     # some proportional area (here 75% or 95%) meet criteria of
    #     # (-9 <= depth <= -3, slope < 40, and habitat is Rock or Coral/Algae).
    #     # 75% is assessed for comparison purposes.
    #     # suitable = read(
    #     #     depth_criteria(src_bathy) .& slope_criteria(src_slope) .& supports_coral(src_benthic)
    #     # )

    #     # See comment above re suitability functions - use of functions breaks `read()`
    #     suitable_flats = read(
    #         (src_geomorphic .∈ [FLAT_IDS]) .&
    #         (src_benthic .∈ [BENTHIC_IDS]) .&
    #         (-9.0 .<= src_bathy .<= -2.0) .&
    #         (0.0 .<= src_slope .<= 40.0)
    #     )

    #     suitable_slopes = read(
    #         (src_geomorphic .∈ [SLOPE_IDS]) .&
    #         (src_benthic .∈ [BENTHIC_IDS]) .&
    #         (-9.0 .<= src_bathy .<= -2.0) .&
    #         (0.0 .<= src_slope .<= 40.0)
    #     )

    #     # Need a copy of raster data type to support writing to `tif`
    #     result_raster = convert.(Int16, copy(suitable_flats))
    #     rebuild(result_raster; missingval=0)

    #     # 85% threshold
    #     # Assess flats
    #     res85 = mapwindow(suitability_func(0.85), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
    #     fpath = joinpath(RESULT_DIR, "$(reg)_suitable_flats_85.tif")
    #     result_raster .= res85
    #     write(fpath, result_raster; force=true)

    #     fpath = joinpath(RESULT_DIR, "$(reg)_grouped_flats_85.tif")
    #     result_raster .= label_components(res85)
    #     write(fpath, result_raster; force=true)

    #     # Assess slopes
    #     res85 = mapwindow(suitability_func(0.85), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
    #     fpath = joinpath(RESULT_DIR, "$(reg)_suitable_slopes_85.tif")
    #     result_raster .= res85
    #     write(fpath, result_raster; force=true)

    #     fpath = joinpath(RESULT_DIR, "$(reg)_grouped_slopes_85.tif")
    #     result_raster .= label_components(res85)
    #     write(fpath, result_raster; force=true)

    #     res85 = nothing
    #     GC.gc()

    #     ####

    #     # 95% threshold
    #     # Assess flats
    #     res95 = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
    #     fpath = joinpath(RESULT_DIR, "$(reg)_suitable_flats_95.tif")
    #     if !isfile(fpath)
    #         result_raster .= res95
    #         write(fpath, result_raster; force=true)
    #     end

    #     fpath = joinpath(RESULT_DIR, "$(reg)_grouped_flats_95.tif")
    #     if !isfile(fpath)
    #         result_raster .= label_components(res95)
    #         write(fpath, result_raster; force=true)
    #     end

    #     # Assess slopes
    #     res95 = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
    #     fpath = joinpath(RESULT_DIR, "$(reg)_suitable_slopes_95.tif")
    #     if !isfile(fpath)
    #         result_raster .= res95
    #         write(fpath, result_raster; force=true)
    #     end
        
    #     fpath = joinpath(RESULT_DIR, "$(reg)_grouped_slopes_95.tif")
    #     if !isfile(fpath)
    #         result_raster .= label_components(res95)
    #         write(fpath, result_raster; force=true)
    #     end

    #     res95 = nothing
    #     GC.gc()

    #     # Apply mode filter (doesn't work as intended)
    #     # mode_res = mapwindow(c_mode, res95, (5, 5))
    #     # result_raster .= label_components(mode_res)
    #     # write(joinpath(RESULT_DIR, "$(reg)_grouped_filtered_suitability.tif"), result_raster; force=true)
    #     # mode_res = nothing

    #     # res75 = nothing
    #     # res95 = nothing
    #     # result_raster = nothing
    # end

    function mask_by_class(raster::AbstractRaster; geo_df::GeoDataFrames.GeoDataFrame, classes::Vector{String})
        # Filter the GeoDataFrame based on the specified classes
        filtered_geo_df = filter(row -> in(row.class, classes), geo_df)
        
        # Mask the raster with the filtered GeoDataFrame
        masked_raster = Rasters.mask(raster; with=filtered_geo_df)
        
        return masked_raster
    end    

    function assess_Allen(LAYERS)         # function adapted from assess_region in 2_apply_analysis.jl
        allen_bathy_path = first(glob("*.tif", joinpath(ALLEN_DATA_DIR, LAYERS[1])))
        allen_bathy = Raster(allen_bathy_path, crs=EPSG(4326), lazy=true)
    
        allen_turbid_path = first(glob("*.tif", joinpath(ALLEN_DATA_DIR, LAYERS[2])))
        allen_turbi = Raster(allen_turbid_path, crs=EPSG(4326), lazy=true)


        allen_benthic_path = first(glob("*.geojson", joinpath(ALLEN_DATA_DIR, "Benthic-Map")))
        allen_benthic = gdf = GDF.read(allen_benthic_path)
    
        allen_geo_path = first(glob("*.geojson", joinpath(ALLEN_DATA_DIR, "Geomorphic-Map")))
        allen_geo = gdf = GDF.read(allen_geo_path)

        # identify flats

        # identify slopes
    
        # # Source image is of 10m^2 pixels
        # # A hectare is 100x100 meters, so we're looking for contiguous areas where
        # # some proportional area (here 75% or 95%) meet criteria of
        # # (-9 <= depth <= -3, slope < 40, and habitat is Rock or Coral/Algae).
        # # 75% is assessed for comparison purposes.
        # # suitable = read(
        # #     depth_criteria(src_bathy) .& slope_criteria(src_slope) .& supports_coral(src_benthic)
        # # )
    
        # See comment above re suitability functions - use of functions breaks `read()`
    
        ## update a/b values ##
        a = -1
        b = 1
        suitable_areas = read(
            (-9.0 .<= allen_bathy .<= -2.0) .&
            (a .<= allen_turbi .<= b) 
            # .& wave data
        )
    
        BENTHICS = ["Rock", "Coral/Algae"]
        FLATS = ["Inner Reef Flat", "Outer Reef Flat", "Plateau"]
        SLOPES = ["Sheltered Reef Slope", "Reef Slope", "Back Reef Slope"]
        ## Geomorphic zones unaccounted for
        # "Deep Lagoon"
        # "Shallow Lagoon"
        # "Reef Crest"
        # "Terrestrial Reef Flat"
        # "Patch Reefs"

        suitable_benthic = mask_by_class(suitable_areas, allen_benthic, BENTHIC_IDS)

        suitable_flats = mask_by_class(suitable_benthic, allen_geo, FLATS)
        # Rasters.mask(suitable_areas; with=[allen_geo.class .== ("inner-reef-flat") .| ("outer-reef-flat") .| ("plateau")])
        
        suitable_slopes = mask_by_class(suitable_benthic, allen_geo, SLOPES)
            # Rasters.mask(suitable_areas; with=[allen_geo.class .== ("sheltered-reef-slope") .| ("reef-slope") .| ("back-reef-slope")])

        ###############################################
        #### RUNS TO HERE ####
        ###############################################

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_flats))
        rebuild(result_raster; missingval=0)
    
        ####
    
        # 95% threshold
        # Assess flats
        res95 = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(ALLEN_RESULT_DIR, "$(reg)_suitable_flats_95.tif")
        if !isfile(fpath)
            result_raster .= res95
            write(fpath, result_raster; force=true)
        end
    
        fpath = joinpath(ALLEN_RESULT_DIR, "$(reg)_grouped_flats_95.tif")
        if !isfile(fpath)
            result_raster .= label_components(res95)
            write(fpath, result_raster; force=true)
        end
    
        # Assess slopes
        res95 = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(ALLEN_RESULT_DIR, "$(reg)_suitable_slopes_95.tif")
        if !isfile(fpath)
            result_raster .= res95
            write(fpath, result_raster; force=true)
        end
        
        fpath = joinpath(ALLEN_RESULT_DIR, "$(reg)_grouped_slopes_95.tif")
        if !isfile(fpath)
            result_raster .= label_components(res95)
            write(fpath, result_raster; force=true)
        end
    
        res95 = nothing
        GC.gc()
    
        ####
        
        # # 85% threshold
        # # Assess flats
        # res85 = mapwindow(suitability_func(0.85), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # fpath = joinpath(RESULT_DIR, "$(reg)_suitable_flats_85.tif")
        # result_raster .= res85
        # write(fpath, result_raster; force=true)
    
        # fpath = joinpath(RESULT_DIR, "$(reg)_grouped_flats_85.tif")
        # result_raster .= label_components(res85)
        # write(fpath, result_raster; force=true)
    
        # # Assess slopes
        # res85 = mapwindow(suitability_func(0.85), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        # fpath = joinpath(RESULT_DIR, "$(reg)_suitable_slopes_85.tif")
        # result_raster .= res85
        # write(fpath, result_raster; force=true)
    
        # fpath = joinpath(RESULT_DIR, "$(reg)_grouped_slopes_85.tif")
        # result_raster .= label_components(res85)
        # write(fpath, result_raster; force=true)
    
        # res85 = nothing
        # GC.gc()
    
        ####
    
        # Apply mode filter (doesn't work as intended)
        # mode_res = mapwindow(c_mode, res95, (5, 5))
        # result_raster .= label_components(mode_res)
        # write(joinpath(RESULT_DIR, "$(reg)_grouped_filtered_suitability.tif"), result_raster; force=true)
        # mode_res = nothing
    
        # res75 = nothing
        # res95 = nothing
        # result_raster = nothing
    end

    # # For some reason, wrapping conditions in functions causes `read()` to fail
    # # Kept here in case it becomes useful later.
    # function depth_criteria(depth)
    #     return (-9 .<= depth .<= -2.0)
    # end
    # function slope_criteria(slope)
    #     return (0.0 .< slope .< 40.0)
    # end
    # function supports_coral(benthic)
    #     # Manually extracted from Raster Attribute Table
    #     # 0x00 = 0 = missing value
    #     # 0x0b = 11 = Sand
    #     # 0x0c = 12 = Rubble
    #     # 0x0d = 13 = Rock
    #     # 0x0f = 15 = Coral/Algae
    #     return (benthic .∈ [[0x0d, 0x0f]])
    # end
    # # Image filtering functions
    # function c_median(x)
    #     if sum(x) == 0
    #         return 0
    #     end

    #     return floor(Int64, median(x))
    # end
    # function c_mode(x)
    #     # Check center of odd-sized square matrix
    #     d = size(x)
    #     centroid = ceil.(Int64, d .* 0.5)
    #     if x[centroid...] == 0
    #         return 0
    #     end

    #     total = sum(x)
    #     if total <= 3
    #         return 0
    #     end

    #     return mode(x)
    # end
end

# @showprogress dt=10 desc="Analyzing..." pmap(assess_region, REGIONS)

@showprogress dt=10 desc="Analyzing..." pmap(assess_Allen, LAYERS)

