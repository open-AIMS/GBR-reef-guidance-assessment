"""
Identify suitable locations based on bathymetric, slope, geomorphic, benthic and wave criteria.
Use MPA raster data to analyze suitable areas.
Output raster files with proportional suitability.
"""

include("common.jl")


"""
    TODO: Replace with polygon-based search
"""
function proportion_suitable(x::BitMatrix, window::Tuple)::Matrix{Int16}
    x′ = zeros(Int16, size(x))

    @floop for row_col in ThreadsX.findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col + window[1], 1)
        x_right = min(col + window[2], size(x, 2))

        y_top = max(row + window[1], 1)
        y_bottom = min(row + window[2], size(x, 1))

        x′[row, col] = Int16(sum(@views x[y_top:y_bottom, x_left:x_right]))
    end

    return x′
end

function load_and_assess(data_path, lb, ub)
    src_data = Raster(data_path; lazy=true)
    in_criteria = (lb .<= src_data .<= ub)

    return in_criteria
end
function load_and_assess(data_path, func)
    src_data = Raster(data_path; lazy=true)

    return func(src_data)
end

function assess_region(reg, port_buffer)
    @info reg

    @info "Reading"
    @time begin
    # Load required prepared raster files for analysis
    bathy_crit = load_and_assess(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"), -9.0, -2.0)
    slope_crit = load_and_assess(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"), 0.0, 40.0)

    benthic_crit = load_and_assess(
        joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
        x -> (x .∈ [MPA_BENTHIC_IDS])
    )

    src_geomorphic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"))
    geomorphic_flat_crit = src_geomorphic .∈ [MPA_FLAT_IDS]
    geomorphic_slope_crit = src_geomorphic .∈ [MPA_SLOPE_IDS]
    src_geomorphic = nothing

    Hs_waves_crit = load_and_assess(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"), 0.0, 1.0)
    Tp_waves_crit = load_and_assess(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"), 0.0, 6.0)
    turbid_crit = load_and_assess(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"), 0.0, 58)
    end

    @info "Assessment"
    @info "Initial assessment"
    @time begin

    # Apply filtering criteria to raster grid
    @time begin
    suitable_areas = (
        bathy_crit .&
        benthic_crit .&
        slope_crit .&
        Hs_waves_crit .&
        Tp_waves_crit .&
        turbid_crit
    )

    if reg == "Townsville-Whitsunday"
        src_rugosity = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"); lazy=true)
        suitable_areas = suitable_areas .& (src_rugosity .<= 6)
        src_rugosity = nothing
    end
    end  # end raster comparison

    # Filter out cells over 200NM from the nearest port
    suitable_areas = filter_distances(suitable_areas, port_buffer)

    # Filter out cells occurring in preservation zones
    GBRMPA_zone_exclusion = GDF.read(joinpath(MPA_OUTPUT_DIR, "GBRMPA_preservation_zone_exclusion.gpkg"))
    region_extent = GI.extent(suitable_areas)
    rst_extent = [
        GI.Polygon(
            create_poly(create_bbox(region_extent.X, region_extent.Y), EPSG(7844))
        )
    ]
    in_region = GO.within.(GBRMPA_zone_exclusion.geometry, rst_extent)
    GBRMPA_zone_exclusion = GBRMPA_zone_exclusion[in_region, :]

    for polygon in GBRMPA_zone_exclusion.geometry
        suitable_areas = mask(suitable_areas; with=polygon, invert=true, boundary=:touches)
    end

    end  # initial assessment

    # Need a copy of raster data type to support writing to `tif`
    # Assess flats
    @info "    Assess flats"
    @time begin

    result_raster = convert.(Int16, copy(suitable_areas))
    rebuild(result_raster; missingval=0)

    suitable_flats = suitable_areas .& geomorphic_flat_crit
    cleaned_flats = remove_orphaned_elements(BitMatrix(suitable_flats.data), 7, (3, 3))
    cleaned_flats = remove_orphaned_elements(cleaned_flats, 70, (9, 9))

    suitable_flats.data .= cleaned_flats
    cleaned_flats = nothing

    # Calculate suitability of 10x10m surroundings of each cell
    res = proportion_suitable(suitable_flats.data, (-4, 5))
    if reg == "Townsville-Whitsunday"
        fpath = joinpath(MPA_ANALYSIS_RESULTS, "$(reg)_suitable_flats_rugosity.tif")
    else
        fpath = joinpath(MPA_ANALYSIS_RESULTS, "$(reg)_suitable_flats.tif")
    end
    result_raster .= res

    end  # end assess flats

    write(fpath, result_raster; force=true)
    suitable_flats = nothing

    # Assess slopes
    @info "    Assess slopes"
    @time begin
    suitable_slopes = suitable_areas .& geomorphic_slope_crit

    # Calculate suitability of 10x10m surroundings of each cell
    res .= proportion_suitable(suitable_slopes.data)
    if reg == "Townsville-Whitsunday"
        fpath = joinpath(MPA_ANALYSIS_RESULTS, "$(reg)_suitable_slopes_rugosity.tif")
    else
        fpath = joinpath(MPA_ANALYSIS_RESULTS, "$(reg)_suitable_slopes.tif")
    end

    end  # assess slopes

    result_raster .= res
    write(fpath, result_raster; force=true)

    suitable_areas = nothing
    suitable_slopes = nothing
    res = nothing

    GC.gc()
end

# Load QLD_ports buffer data
port_buffer = GDF.read(joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg"))[:,:geometry]

@showprogress dt = 10 desc = "Analyzing..." map(x -> assess_region(x, port_buffer), REGIONS)
