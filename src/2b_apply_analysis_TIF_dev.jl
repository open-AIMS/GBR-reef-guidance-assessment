"""Identify suitable locations per region."""

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

function mask_by_class(raster::AbstractRaster; geo_df::Any, classes::Vector{String})
    # Filter  GDF based on the specified classes
    filtered_geo_df = filter(row -> in(row.class, classes), geo_df)
    
    # Mask raster with filtered GeoDataFrame
    masked_raster = Rasters.mask(raster; with=filtered_geo_df)    
    return masked_raster
end    

function assess_allen()         # function adapted from assess_region in 2_apply_analysis.jl
    allen_bathy_path = first(glob("*.tif", joinpath(ALLEN_DATA_DIR, "Bathymetry---composite-depth")))
    allen_bathy = Raster(allen_bathy_path, crs=EPSG(4326), lazy=true)

    allen_turbid_path = first(glob("*.tif", joinpath(ALLEN_DATA_DIR, "Turbidity-Q3-2023")))
    allen_turbi = Raster(allen_turbid_path, crs=EPSG(4326), lazy=true)

    ## update a/b values ##
    suitable_areas = read(
        (-9.0 .<= allen_bathy .<= -2.0) .&  # depth criteria
        (allen_turbi .<= 52)    # LOW Turbidity = “52” in the turbidity maps = “5.2 FNU”
        # .& wave data
    )    
    allen_bathy = nothing
    GC.gc()
    allen_turbi = nothing
    GC.gc()

    allen_benthic_path = first(glob("*.geojson", joinpath(ALLEN_DATA_DIR, "Benthic-Map")))
    allen_benthic = GDF.read(allen_benthic_path)
    suitable_benthic = mask_by_class(suitable_areas, allen_benthic, BENTHIC_ALLEN)
    allen_benthic = nothing
    GC.gc()

    allen_geo_path = first(glob("*.geojson", joinpath(ALLEN_DATA_DIR, "Geomorphic-Map")))
    allen_geo = GDF.read(allen_geo_path)
    suitable_flats = mask_by_class(suitable_benthic, allen_geo, FLATS_ALLEN)        
    suitable_slopes = mask_by_class(suitable_benthic, allen_geo, SLOPES_ALLEN)
    allen_geo = nothing
    GC.gc()

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

assess_allen()

