using Rasters
import GeoDataFrames as GDF
import ArchGDAL as AG
using DataFrames
using ProgressMeter
using Glob
using Makie


const DATA_DIR = "C:/Users/rlippman/Documents/development/ADRIA_data/spatial_datasets/GBR-Bathy10m"
const ALLEN_DIR = joinpath(DATA_DIR, "..", "AllenAtlas_GBR-20231118074407")

const criteria =  ["Bathymetry---composite-depth", 
                    "Turbidity-Q3-2023"]          #readdir(DATA_DIR)    # Criteria = all files/folders in the directory


allen_bathy_path = "$(ALLEN_DIR)/$(criteria[1])/bathymetry_0.tif"
allen_bathy = Raster(allen_bathy_path, crs=EPSG(4326), lazy=true)

allen_turbid_path = "$(ALLEN_DIR)/$(criteria[2])/turbidity-quarterly_0.tif"
allen_turbi = Raster(allen_turbid_path, crs=EPSG(4326), lazy=true)

const RESULT_DIR = joinpath(ALLEN_DIR, "outputs")

# # Is this equivalent to Allen_Atlas Reef-Extent??
# Get polygon of management areas
# region_path = joinpath(
#     DATA_DIR,
#     "zones",
#     "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
# )
# region_features = GDF.read(region_path)

function assess_Allen(criteria)         # function adapted from assess_region in 2_apply_analysis.jl
    allen_bathy_path = first(glob("*.tif", joinpath(ALLEN_DIR, criteria[1])))
    allen_bathy = Raster(allen_bathy_path, crs=EPSG(4326), lazy=true)

    allen_turbid_path = first(glob("*.tif", joinpath(ALLEN_DIR, criteria[2])))
    allen_turbi = Raster(allen_turbid_path, crs=EPSG(4326), lazy=true)

    # # Load pre-prepared data - pre-prepared in 1b_prep_data.jl
    # ##allen_benthic_path = "../figs/$(reg)_benthic.tif"
    # allen_benthic = Raster(allen_benthic_path, lazy=true)
    # ##allen_geomorphic_path = "../figs/$(reg)_geomorphic.tif"
    # allen_geo = Raster(allen_geomorphic_path, lazy=true)


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
        (allen_bathy) .&
        (allen_turbi) .&
        (allen_geo .∈ [FLAT_IDS]) .&            # Inner Reef Flat, Outer Reef Flat, Plateau
        (allen_benthic .∈ [BENTHIC_IDS])        # Rock, Coral/Algae
    )

    suitable_slopes = read(
        (allen_bathy) .&
        (allen_turbi) .&
        (allen_geo .∈ [SLOPE_IDS]) .&           # Sheltered Reef Slope, Reef Slope, Back Reef Slope
        (allen_benthic .∈ [BENTHIC_IDS])        # Rock, Coral/Algae
    )

    # Need a copy of raster data type to support writing to `tif`
    result_raster = convert.(Int16, copy(suitable_flats))
    rebuild(result_raster; missingval=0)

    ####

    # 95% threshold
    # Assess flats
    res95 = mapwindow(suitability_func(0.95), suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
    fpath = joinpath(RESULT_DIR, "$(reg)_suitable_flats_95.tif")
    if !isfile(fpath)
        result_raster .= res95
        write(fpath, result_raster; force=true)
    end

    fpath = joinpath(RESULT_DIR, "$(reg)_grouped_flats_95.tif")
    if !isfile(fpath)
        result_raster .= label_components(res95)
        write(fpath, result_raster; force=true)
    end

    # Assess slopes
    res95 = mapwindow(suitability_func(0.95), suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
    fpath = joinpath(RESULT_DIR, "$(reg)_suitable_slopes_95.tif")
    if !isfile(fpath)
        result_raster .= res95
        write(fpath, result_raster; force=true)
    end
    
    fpath = joinpath(RESULT_DIR, "$(reg)_grouped_slopes_95.tif")
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


# output_csv_path_bathy = joinpath(RESULT_DIR, "allen_bathy.csv")
# output_csv_path_turbid = joinpath(RESULT_DIR, "allen_turbid.csv")


# Define function to export raster to CSV
function export_raster_to_csv(raster, output_csv_path)
    nrows, ncols = size(raster)
    println(size(raster))
    println("Size of raster: ($(ncols), $(nrows))")
    chunk_size = 10  # Adjust the chunk size as needed
    nrows = 1000
    ncols = 10
    println("Size of considered raster: ($(ncols), $(nrows))")
    # println(view(raster, 1:5, :))
    for i in 1:chunk_size:nrows       # Loop through rows
        nrows_chunk = min(chunk_size, nrows - i + 1)
        chunk = view(raster, i:i+nrows_chunk-1, :)
        
        df_chunk = DataFrame(chunk)
        
        if i == 1
            CSV.write(output_csv_path, df_chunk)
        else
            CSV.write(output_csv_path, df_chunk; append=true, header=false)
        end
    end
end
export_raster_to_csv(allen_bathy, output_csv_path_bathy)
export_raster_to_csv(allen_turbid, output_csv_path_turbid)
