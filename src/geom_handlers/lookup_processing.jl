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

"""
    valid_lookup(raster_files::NamedTuple, valid_areas_file::String, dst_file::String)::Nothing

Create a lookup table of valid data pixels for fast querying of data layers.

# Arguments
- `raster_files` : NamedTuple containing the file path for each criteria raster file.
- `valid_areas_file` : Path for file containing target valid areas (slopes or flats).
- `dst_file` : Path to write parquet lookup file to.
"""
function valid_lookup(raster_files::NamedTuple, valid_areas_file::String, dst_file::String)::Nothing
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    # Create stack of prepared data
    rst_stack = RasterStack(raster_files; lazy=true)

    # Create lookup of valid data
    valid_areas = Raster(valid_areas_file)
    _valid = sparse(boolmask(valid_areas).data)
    valid_areas = nothing
    force_gc_cleanup()

    col_names = vcat(:geometry, :lon_idx, :lat_idx, keys(raster_files)...)
    area_values = stack_values(_valid, rst_stack)
    area_store = geoparquet_df!(area_values, col_names)
    GP.write(dst_file, area_store, (:geometry, ))

    area_store = nothing
    area_values = nothing
    _valid = nothing
    rst_stack = nothing
    force_gc_cleanup(; wait_time=4)

    return nothing
end
