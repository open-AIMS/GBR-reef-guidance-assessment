"""
Identify suitable locations based on bathymetric, slope, geomorphic, benthic and wave criteria.
Use MPA raster data to analyze suitable areas.
Output raster files with proportional suitability.
"""

include("common.jl")


"""
    TODO: Replace with polygon-based search
"""
function proportion_suitable(x::BitMatrix)::Matrix{Int16}
    x′ = zeros(Int16, size(x))

    @floop for row_col in ThreadsX.findall(x)
        (row, col) = Tuple(row_col)
        x_left = max(col - 4, 1)
        x_right = min(col + 4, size(x, 2))

        y_top = max(row - 4, 1)
        y_bottom = min(row + 4, size(x, 1))

        x′[row, col] = Int16(sum(@views x[y_top:y_bottom, x_left:x_right]))
    end

    return x′
end

"""
    port_buffer_mask(gdf::DataFrame, dist::Float64; unit::String="NM")

Create a masking buffer around indicated port locations.

# Arguments
- `gdf` : GeoDataFrame of port locations (given as long/lat points)
- `dist` : distance from port in degrees (deg), kilometers (km), or nautical miles (NM; default)
- `unit` : unit `dist` is in
"""
function port_buffer_mask(gdf::DataFrame, dist::Float64; unit::String="NM")
    # Determine conversion factor (nautical miles or kilometers)
    conv_factor = 1.0
    if unit == "NM"
        conv_factor = 60.0  # 60 NM = 1 degree
    elseif unit == "km"
        conv_factor = 111.0  # 111 km = 1 degree
    elseif unit != "deg"
        error("Unknown distance unit requested. Can only be one of `NM` or `km` or `deg`")
    end

    ports = gdf.geometry  # TODO: Replace with `GI.geometrycolumns()`

    # Make buffer around ports
    buffered_ports = GO.buffer.(ports, dist / conv_factor)

    # Combine all geoms into one
    port_mask = reduce((x1, x2) -> LibGEOS.union(x1, x2), buffered_ports)

    return port_mask
end

"""
    filter_distances(
        target_rast::Raster,
        dist_buffer
    )::Raster

Apply a mask to exclude pixels that are outside the indicated distance buffer(s).

`target_rast` and the `dist_buffer` should be in the same CRS (e.g., EPSG:7844 / GDA2020).

# Arguments
- `target_rast` : Raster of suitable pixels (Bool) to filter pixels from.
- `dist_buffer` : Buffer geometry to use as the mask.

# Returns
- `within_dist` : Masked raster of indicating pixels that are within the target distance.
"""
function filter_distances(target_rast::Raster, dist_buffer)::Raster
    # Mask out areas outside considered distance from port
    return mask(Raster(target_rast; missingval=0); with=dist_buffer)
end

function load_and_assess(data_path, lb, ub)
    src_data = Raster(data_path)
    in_criteria = (lb .<= src_data .<= ub)

    return in_criteria
end
function load_and_assess(data_path, func)
    src_data = Raster(data_path)

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

    # Filter out cells over 200NM from the nearest port
    bathy_crit = filter_distances(bathy_crit, port_buffer)

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
        src_rugosity = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"))
        suitable_areas .= suitable_areas .& (src_rugosity .< 6)
        src_rugosity = nothing
    end
    end  # end raster comparison

    end  # initial assessment

    # Need a copy of raster data type to support writing to `tif`
    # Assess flats
    @info "    Assess flats"
    @time begin

    result_raster = convert.(Int16, copy(suitable_areas))
    rebuild(result_raster; missingval=0)

    suitable_flats = suitable_areas .& geomorphic_flat_crit

    # Calculate suitability of 10x10m surroundings of each cell
    res = proportion_suitable(suitable_flats.data)
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

# Load QLD_ports data
port_locs = GDF.read("$(PORT_DATA_DIR)/ports_QLD_merc.shp")
port_locs.geometry = AG.reproject(
    port_locs.geometry,
    crs(port_locs[1, :geometry]),
    GDA2020_crs;
    order=:trad
)

port_buffer = port_buffer_mask(port_locs, 200.0)

@showprogress dt = 10 desc = "Analyzing..." map(x -> assess_region(x, port_buffer), REGIONS)
