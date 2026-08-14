"""
    stack_values(valid_mask, rst_stack; band_height::Int=2048)

Extract values at specific lon/lat coordinates from a raster stack.

Reads each layer in GDAL-aligned row bands (height a multiple of 256, matching the
256x256 COG block layout of the criteria rasters) instead of materialising each
layer's full bounding box, and scatters values directly into preallocated typed
columns instead of a boxed `Matrix{Any}`.

# Notes
Currently expects the raster to have the default X/Y dimensions set.

# Arguments
- `valid_mask` : mask (anything `SparseArrays.SparseMatrixCSC` can be built from)
  indicating locations of valid data
- `rst_stack` : raster stack to extract data from
- `band_height` : number of Y rows to read per band; must be a multiple of 256 to
  stay aligned with the criteria layers' COG block size

# Returns
Tuple of `(value_cols, lon_idx_col, lat_idx_col, centroid_lons, centroid_lats)`.
`value_cols` is a `NamedTuple` of one typed vector per stack layer (eltype matches
`eltype(rst_stack[name])`), `lon_idx_col`/`lat_idx_col` are `Int64` source-pixel
indices, and `centroid_lons`/`centroid_lats` are centroid-adjusted coordinates - all
five are in the same row order: sorted `(x, y)` lexicographically (X-major), matching
the pipeline's existing Arrow output. Callers that also need lon/lat columns (e.g.
`valid_lookup`) should reuse these rather than re-deriving them.
"""
function stack_values(valid_mask, rst_stack; band_height::Int=2048)
    @assert band_height % 256 == 0 "band_height must be a multiple of 256 to stay aligned with COG blocks"

    # Calculate cell sizes
    x_res = abs(step(lookup(rst_stack, X)))
    y_res = abs(step(lookup(rst_stack, Y)))

    # Collect locations in lat/longs
    lons = collect(lookup(rst_stack, X))
    lats = collect(lookup(rst_stack, Y))

    # Walk the mask's nonzeros in native column-major (Y-major) order: for a matrix
    # with dims (X, Y), that is every valid pixel in image row (Y) 1, then row 2, and
    # so on - the same order the GDAL row-band reads below visit pixels in, so no
    # separate re-ordering pass is needed to line the two up.
    csc_mask = SparseArrays.SparseMatrixCSC(valid_mask)
    x_idx, y_idx, _ = findnz(csc_mask)
    n_valid = length(x_idx)

    # Output row order must stay X-major (sorted (x, y) lexicographically) to match
    # the pipeline's existing Arrow output (every valid pixel in column 1, then
    # column 2, ...). `perm[k]` is the destination row for the k-th pixel visited
    # during the Y-major scan below.
    #
    # `sortperm` alone gives *source* indices (coord_list[p[r]] is the pixel of rank
    # r), not destinations. Scattering with `dest[perm[k]] = src[k]` needs each
    # pixel's own rank, which is `invperm(sortperm(coord_list))`.
    perm = invperm(sortperm(collect(zip(x_idx, y_idx))))

    centroid_lons = Vector{Float64}(undef, n_valid)
    centroid_lats = Vector{Float64}(undef, n_valid)
    lon_idx_col = Vector{Int64}(undef, n_valid)
    lat_idx_col = Vector{Int64}(undef, n_valid)
    for k in 1:n_valid
        dest = perm[k]
        centroid_lons[dest] = lons[x_idx[k]] + (x_res / 2)
        centroid_lats[dest] = lats[y_idx[k]] + (y_res / 2)
        lon_idx_col[dest] = x_idx[k]
        lat_idx_col[dest] = y_idx[k]
    end

    # Typed column store: one preallocated vector per layer, eltype derived from the
    # layer itself (Benthic/Geomorphic are Int8 class IDs, the rest are Float32).
    stack_names = names(rst_stack)
    value_cols = NamedTuple{Tuple(stack_names)}(
        Tuple(Vector{eltype(rst_stack[name])}(undef, n_valid) for name in stack_names)
    )

    ny = length(lats)
    k = 1
    y_start = 1
    while y_start <= ny
        y_end = min(y_start + band_height - 1, ny)

        # y_idx is non-decreasing (CSC column-major order), so the pixels belonging
        # to this band form one contiguous run starting at k.
        band_start = k
        while k <= n_valid && y_idx[k] <= y_end
            k += 1
        end
        band_end = k - 1

        if band_start <= band_end
            for stack_name in stack_names
                # Plain getindex on the lazy raster is a windowed GDAL read; `[:, :]`
                # forces that window to be read into memory.
                band_data = rst_stack[stack_name][:, y_start:y_end][:, :]
                col = value_cols[stack_name]
                for j in band_start:band_end
                    col[perm[j]] = band_data[x_idx[j], y_idx[j]-y_start+1]
                end
            end
        end

        y_start = y_end + 1
        force_gc_cleanup()
    end

    return value_cols, lon_idx_col, lat_idx_col, centroid_lons, centroid_lats
end

"""
    valid_lookup(raster_files::NamedTuple, valid_areas_file::String, dst_file::String)::Nothing

Create a lookup table of valid data pixels for fast querying of data layers.

# Arguments
- `raster_files` : NamedTuple containing the file path for each criteria raster file.
- `valid_areas_file` : Path for file containing target valid areas (slopes or flats).
- `dst_file` : Path to write Arrow lookup file to.
"""
function valid_lookup(raster_files::NamedTuple, valid_areas_file::String, dst_file::String)::Nothing
    lookup_sources = (collect(values(raster_files))..., valid_areas_file)
    return skip_if_exists(dst_file; sources=lookup_sources) do
        # Create stack of prepared data
        rst_stack = RasterStack(raster_files; lazy=true)

        # Create lookup of valid data
        valid_areas = Raster(valid_areas_file; lazy=true)
        _valid = ExtendableSparseMatrix(sparse(boolmask(valid_areas).data))
        valid_areas = nothing
        force_gc_cleanup()

        value_cols, lon_idx_col, lat_idx_col, centroid_lons, centroid_lats =
            stack_values(_valid, rst_stack)

        area_store = DataFrame(; lon_idx=lon_idx_col, lat_idx=lat_idx_col, value_cols...)

        # Reuse the centroid coordinates already computed in `stack_values` instead of
        # re-deriving them from a geometry column via a redundant geometry decode pass.
        area_store[!, :lons] = centroid_lons
        area_store[!, :lats] = centroid_lats

        Arrow.write(dst_file, area_store; compress=:zstd)

        # Write a JSON sidecar of per-criterion (min, max) bounds alongside the Arrow
        # lookup, for lazy-scoped readers that need to know a column's value range
        # without loading the full Arrow table.
        bounds_path = replace(
            dst_file, "_valid_slopes_lookup.arrow" => "_valid_slopes_bounds.json"
        )
        skip_cols = Set([:lon_idx, :lat_idx, :lons, :lats])
        skip_if_exists(bounds_path; label="Bounds sidecar") do
            bounds_dict = Dict{String,Any}()
            for col in names(area_store)
                sym = Symbol(col)
                sym in skip_cols && continue
                vals = area_store[!, col]
                if nonmissingtype(eltype(vals)) <: Real
                    non_missing = skipmissing(vals)
                    if !isempty(non_missing)
                        bounds_dict[col] = Dict(
                            "min" => Float64(minimum(non_missing)),
                            "max" => Float64(maximum(non_missing))
                        )
                    end
                end
            end
            open(bounds_path, "w") do io
                JSON3.write(io, bounds_dict)
            end

            return nothing
        end

        area_store = nothing
        value_cols = nothing
        _valid = nothing
        rst_stack = nothing
        force_gc_cleanup()

        return nothing
    end
end
