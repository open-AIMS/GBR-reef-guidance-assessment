"""
Prepare data for analysis by processing MPA files for each GBRMPA management region.

Crop GBR-wide GBRMPA rasters into management regions.
Reproject all data from WGS84 / UTM Zone 54 - 56 into consistent CRS (GDA2020).
Ensure all rasters are the same size/shape for each region of interest with the same
values used to indicate no data.

The general approach is to crop down to an extent of the region, then trim/mask to just the
areas of interest. Processing data in this way minimized the amount of data being handled.

Still, the largest chunk of time is spent resampling to the same projections and writing
data to file (as it takes time to compress the data). The datasets can be in the GBs without
compression.
"""

include("common.jl")
include("raster_processing.jl")

# This processing step requires all the memory available so we actually remove the workers
# that were just set up.
# This can be commented out if working on a system with large enough amounts of memory
rmprocs(workers()[2:end]...)
GC.gc()

using SparseArrays

# Loading regions_4326 for cropping of vector and raster data.
regions_4326 = GDF.read(REGION_PATH_4326)

# 1. Process GBRMPA zoning geopackage to select only zones for site exclusion
if !isfile(joinpath(MPA_OUTPUT_DIR, "GBRMPA_zone_exclusion.gpkg"))
    GBRMPA_zoning_poly = GDF.read(
        joinpath(
            GDA2020_DATA_DIR,
            "Great_Barrier_Reef_Marine_Park_Zoning_20_4418126048110066699.gpkg"
        )
    )
    GBRMPA_zoning_poly = GBRMPA_zoning_poly[GBRMPA_zoning_poly.TYPE .∈ [MPA_EXCLUSION_ZONES], :]
    rename!(GBRMPA_zoning_poly, :SHAPE => :geometry)

    GDF.write(
        joinpath(MPA_OUTPUT_DIR, "GBRMPA_zone_exclusion.gpkg"),
        GBRMPA_zoning_poly;
        crs=EPSG_7844
    )
end

if !isfile(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))
    port_locs = GDF.read("$(PORT_DATA_DIR)/ports_QLD_merc.shp")
    port_locs.geometry = AG.reproject(
        port_locs.geometry,
        crs(port_locs[1, :geometry]),
        EPSG_7844;
        order=:trad
    )

    GDF.write(
        joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"),
        port_locs;
        crs=EPSG_7844
    )
end

if !isfile(joinpath(MPA_OUTPUT_DIR, "ports_buffer.gpkg"))
    port_locs = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))

    port_buffer = port_buffer_mask(port_locs, 200.0, unit="NM")
    port_buffer = DataFrame(; Name="ports_buffer", geometry=port_buffer)
    GDF.write(
        joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg"),
        port_buffer;
        crs=EPSG_7844
    )
    port_buffer = nothing
end


# 2. Process MPA files to represent GBRMPA regions in GDA2020 projection

"""
    stack_values(valid_mask, rst_stack)

Extract values at specific lon/lat coordinates from a raster stack.
Loads each stack into memory and extracts values a layer at a time.

# Notes
Currently expects the raster to have the default X/Y dimensions set.

# Arguments
- `valid_mask` : mask indicating locations of valid data
- `rst_stack` : raster stack to extract data from

# Returns
Tables.jl-compatible vector of named tuples (to build a dataframe with)
"""
function stack_values(valid_mask, rst_stack)
    # Collect locations in lat/longs
    lons = collect(lookup(rst_stack, X))
    lats = collect(lookup(rst_stack, Y))

    sorted_valid_idx = sort(Tuple.(findall(valid_mask)))
    lon_lats = collect(zip(lons[first.(sorted_valid_idx)], lats[last.(sorted_valid_idx)]))

    # Create store, with three additional columns to make space for geometry, lon/lat index
    v_store = Matrix(undef, length(lon_lats), length(names(rst_stack))+3)
    for (i, stack_name) in enumerate(names(rst_stack))
        # Read in valid subset
        rst_tmp = read(
            view(
                rst_stack[stack_name],
                sort(unique(first.(sorted_valid_idx))),
                sort(unique(last.(sorted_valid_idx)))
            )
        )

        get_index = i == 1
        extracted = extract(rst_tmp, lon_lats; index=get_index)
        if get_index
            # Returned indices are relative to the view, not the source raster
            # So we jump through some hoops to obtain the canonical indices.
            inds = getfield.(extracted, :index)
            true_lon_inds = lookup(rst_tmp, X).data.indices[1][first.(inds)]
            true_lat_inds = lookup(rst_tmp, Y).data.indices[1][last.(inds)]

            v_store[:, 1] .= getfield.(extracted, :geometry)
            v_store[:, 2] .= true_lon_inds
            v_store[:, 3] .= true_lat_inds
            v_store[:, 4] .= getfield.(extracted, stack_name)
        else
            v_store[:, i+3] .= getfield.(extracted, stack_name)
        end

        rst_tmp = nothing
        extracted = nothing
        force_gc_cleanup(; wait_time=5)
    end

    return v_store
end

"""
    geoparquet_df!(store_values::Matrix, col_names::Vector{Symbol})::DataFrame

Create a GeoParquet-compatible dataframe by assigning correct type information
for each column.

# Arguments
- `store_values` : Values to put into store
- `col_names` : column names to use
"""
function geoparquet_df!(store_values::Matrix, col_names::Vector{Symbol})::DataFrame
    store = try
        DataFrame(store_values, col_names)
    catch
        @warn "Assuming values are in compatible namedtuple"
        DataFrame(store_values)
    end

    for (i, col) in enumerate(eachcol(store))
        store[!, i] = convert.(typeof(store[1, i]), col)
    end

    return store
end

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    # Process MPA bathymetry and slope raster files
    raw_bathy_fn = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    target_bathy_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif")
    process_MPA_bottom_raster(raw_bathy_fn, target_bathy_fn, bathy_gda2020, -9999.0)

    raw_slope_fn = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
    target_slope_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif")
    process_MPA_bottom_raster(raw_slope_fn, target_slope_fn, bathy_gda2020, -9999.0)

    # Process GBR-wide and rugosity raster data
    # Load bathymetry data to provide corresponding spatial extent
    bathy_gda2020 = Raster(base_bathy_fn; lazy=true)

    raw_benthic_fn = "$(MPA_DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
    target_benthic_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif")
    target_benthic = trim_extent_region(
        raw_benthic_fn,
        EPSG_4326,
        regions_4326[reg_idx_4326, :geometry]
    )
    resample_and_write(target_benthic, bathy_gda2020, target_benthic_fn)
    target_benthic = nothing
    force_gc_cleanup()

    raw_geomorphic_fn = "$(MPA_DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
    target_geomorphic_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif")
    target_geomorphic = trim_extent_region(
        raw_geomorphic_fn,
        EPSG_4326,
        regions_4326[reg_idx_4326, :geometry]
    )
    resample_and_write(target_geomorphic, bathy_gda2020, target_geomorphic_fn)
    target_geomorphic = nothing
    force_gc_cleanup()

    raw_turbid_fn = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
    target_turbid_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif")
    target_turbid = trim_extent_region(
        raw_turbid_fn,
        EPSG_4326,
        regions_4326[reg_idx_4326, :geometry]
    )
    resample_and_write(target_turbid, bathy_gda2020, target_turbid_fn)
    target_turbid = nothing
    force_gc_cleanup()

    if reg == "Townsville-Whitsunday"
        raw_rugosity_fn = joinpath(RUG_DATA_DIR, "std25_Rugosity_Townsville-Whitsunday.tif")
        target_rugosity_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif")
        target_rugosity = trim_extent_region(
            raw_rugosity_fn,
            EPSG_4326,
            regions_4326[reg_idx_4326, :geometry]
        )
        resample_and_write(target_rugosity, bathy_gda2020, target_rugosity_fn)
        target_rugosity = nothing
        force_gc_cleanup()
    end

    # Process wave raster data
    # Use bathy dataset as a template for writing netCDF data to geotiff
    src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    rst_template = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))
    target_waves_Hs_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif")
    process_wave_data(waves_Hs_path, target_waves_Hs_fn, :Hs90, rst_template, bathy_gda2020, -9999.0)

    waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Tp", reg)))
    target_waves_Tp_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif")
    process_wave_data(waves_Tp_path, target_waves_Tp_fn, :Tp90, rst_template, bathy_gda2020, -9999.0)

    # Find locations containing valid data
    criteria_paths = Dict(
        "bathy_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"),
        "slope_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"),
        "benthic_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
        "geomorph_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"),
        "waves_Hs_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"),
        "waves_Tp_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"),
        "turbid_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif")
    )
    if reg == "Townsville-Whitsunday"
        push!(criteria_paths, "rugosity_fn" => joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"))
    end

    valid_slopes_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif")
    find_valid_locs(
        criteria_paths,
        MPA_BENTHIC_IDS,
        MPA_SLOPE_IDS,
        7, (3,3), 70, (9,9),
        valid_slopes_fn,
        reg
    )
    valid_flats_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif")
    find_valid_locs(
        criteria_paths,
        MPA_BENTHIC_IDS,
        MPA_FLAT_IDS,
        7, (3,3), 70, (9,9),
        valid_flats_fn,
        reg
    )

    # Calculate distance to nearest port
    port_buffer = GDF.read(joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg"))
    port_points = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))

    port_dist_slopes_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_slopes.tif")
    distance_raster(
        valid_slopes_fn,
        port_buffer,
        port_points,
        -9999.0,
        port_dist_slopes_fn,
        "NM"
    )

    port_dist_flats_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_flats.tif")
    distance_raster(
        valid_flats_fn,
        port_buffer,
        port_points,
        -9999.0,
        port_dist_flats_fn,
        "NM"
    )

    # Create lookup tables to support fast querying
    if !isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes_lookup.parq"))
        # Create stack of prepared data
        # TODO: These paths should be generated elsewhere...
        raster_files = (
            Depth=joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"),
            Slope=joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"),
            Benthic=joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
            Geomorphic=joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"),
            WavesHs=joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"),
            WavesTp=joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"),
            Turbidity=joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"),
            PortDistSlopes=port_dist_slopes_fn,
            PortDistFlats=port_dist_flats_fn
        )

        rst_stack = RasterStack(raster_files; lazy=true)

        # Create lookup of valid slope data
        valid_slopes = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif"))
        _valid = sparse(boolmask(valid_slopes).data)
        valid_slopes = nothing
        force_gc_cleanup()

        col_names = vcat(:geometry, :lon_idx, :lat_idx, keys(raster_files)...)
        slope_values = stack_values(_valid, rst_stack)
        slope_store = geoparquet_df!(slope_values, col_names)
        GP.write(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes_lookup.parq"), slope_store, (:geometry, ))

        slope_store = nothing
        slope_values = nothing
        force_gc_cleanup()

        # Create lookup of valid flat data
        valid_flats = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif"))
        _valid = sparse(boolmask(valid_flats).data)
        valid_flats = nothing
        force_gc_cleanup()

        flat_values = stack_values(_valid, rst_stack)
        flat_store = geoparquet_df!(flat_values, col_names)
        GP.write(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats_lookup.parq"), flat_store, (:geometry, ))

        valid_flats = nothing
        flat_store = nothing
        flat_values = nothing
        _valid = nothing
        rst_stack = nothing
        force_gc_cleanup(; wait_time=4)
    end
end
