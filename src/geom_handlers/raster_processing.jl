using
    Rasters,
    ImageMorphology,
    DataFrames,
    NearestNeighbors

import ArchGDAL as AG
import GeoInterface as GI
import GeoFormatTypes as GFT
import NCDatasets


"""
    write_cog(file_path::String, data::Raster, config::Dict)::Nothing

Write out a COG using common options.

# Arguments
- `file_path` : Path to write data out to
- `data` : Raster data to write out
- `tile_size` : Size of tiles to use in the COG (x,y dimension)
- `num_threads` : Number of threads to use for writing
"""
function write_cog(
    file_path::String, data::Raster; tile_size::Tuple{Integer}=(256,), num_threads::Integer=4
)::Nothing
    Rasters.write(
        file_path,
        data;
        ext=".tiff",
        source="gdal",
        driver="COG",
        options=Dict{String,String}(
            "COMPRESS" => "DEFLATE",
            "SPARSE_OK" => "TRUE",
            "OVERVIEW_COUNT" => "5",
            "BLOCKSIZE" => string(first(tile_size)),
            "NUM_THREADS" => string(num_threads)
        ),
        force=true
    )

    return nothing
end

"""
    resample_to_disk(
        input_raster,
        dst_file::String;
        to=nothing,
        crs=nothing,
        method::Symbol
    )::Nothing

Resample `input_raster` directly to disk (low memory use, but produces an
uncompressed intermediate file), then recompress it into a COG at `dst_file`.

Reading the intermediate file lazily and writing it out with `write_cog`
keeps memory use bounded (GDAL streams block-by-block) while still producing
a compressed result - replicating what the (now removed) ConnectedSystems
fork of Rasters.jl did by default for disk-based `resample()` calls
(see https://github.com/rafaqz/Rasters.jl/issues/706).

# Arguments
- `input_raster` : Raster to resample.
- `dst_file` : Path to write the compressed, resampled result to.
- `to` : Template raster to resample to (mutually exclusive with `crs`).
- `crs` : Target CRS to resample to (mutually exclusive with `to`).
- `method` : Resampling interpolation method.
"""
function resample_to_disk(
    input_raster,
    dst_file::String;
    to=nothing,
    crs=nothing,
    method::Symbol
)::Nothing
    tmp_file = tempname() * ".tif"
    try
        if !isnothing(to)
            resample(input_raster; to=to, method=method, filename=tmp_file)
        else
            resample(input_raster; crs=crs, method=method, filename=tmp_file)
        end
        write_cog(dst_file, Raster(tmp_file; lazy=true))
    finally
        isfile(tmp_file) && rm(tmp_file)
    end

    return nothing
end

"""
    set_consistent_missingval!(raster, val)

Replace value used to indicate no data, and return a Raster type with this value set.

## Note
`raster` should cover the smallest possible extent to reduce processing time.

# Arguments
- `raster` : data to reconfigure
- `val` : Value to use to indicate no data

# Returns
Raster
"""
function set_consistent_missingval!(raster, val)
    # Use non-mutating replace_missing: for disk-backed rasters (lazy=true), Rasters.jl
    # keeps this lazy internally, whereas replace_missing! does an in-place `A .= ...`
    # which requires a writable, fully-realised array and isn't safe/valid for a
    # read-only disk-backed source.
    #
    # Coerce `val` to the raster's own element type first. If `val`'s type differs from
    # the raster's (e.g. a Float64 literal against Float32 data), replace_missing's lazy
    # branch for disk-backed rasters becomes type-unstable internally and DiskArrays
    # widens the resulting eltype to an abstract type (e.g. AbstractFloat) rather than a
    # concrete one - which then breaks downstream GDAL writes/resampling with
    # `convert(GDALDataType, AbstractFloat)` has no method. Matching the type keeps the
    # eltype concrete.
    val = convert(nonmissingtype(eltype(raster)), val)
    raster = replace_missing(raster, val)
    return Raster(raster; missingval=val)
end

"""
    remove_orphaned_elements(rst_mask::BitMatrix, min_cluster_size::Int, box_size::Tuple{Int64,Int64})

Cleans up valid pixels that are by themselves and not worth including in later assessments.

Runs union-find over the set (valid) pixels only, rather than `ImageMorphology.label_components`'s
dense label array - `rst_mask` is typically ~0.5% dense, so this avoids allocating and scanning a
dense `Array{Int}` covering every pixel in the raster (background included) just to label a sparse
handful of foreground components.

# Arguments
- `rst_mask` : Mask of valid raster locations
- `min_cluster_size` : Number of elements that need to be clustered together to be kept
- `box_size` : area to search around center pixel (width, height). Must be odd numbers.
"""
function remove_orphaned_elements(rst_mask::BitMatrix, min_cluster_size::Int, box_size::Tuple{Int64,Int64})
    nx, ny = size(rst_mask)
    coords = findall(rst_mask)
    n = length(coords)

    # Union-find over set-pixel indices (1:n), keyed by linear pixel position so
    # neighbours can be looked up in O(1) without a dense grid.
    pos_to_idx = Dict{Int,Int}()
    sizehint!(pos_to_idx, n)
    for (i, c) in enumerate(coords)
        pos_to_idx[c[1]+(c[2]-1)*nx] = i
    end

    parent = collect(1:n)
    function find_root(a::Int)
        while parent[a] != a
            parent[a] = parent[parent[a]]
            a = parent[a]
        end
        return a
    end

    # `box_size` is the same structuring element `strel_box` builds: two pixels are
    # in the same component if one lies anywhere within the other's `box_size` box,
    # not just its immediate 8-neighbourhood (this matters for the (9,9) second pass).
    rx = box_size[1] ÷ 2
    ry = box_size[2] ÷ 2
    for (i, c) in enumerate(coords)
        x, y = c[1], c[2]
        for dy in 0:ry, dx in (dy == 0 ? (0:rx) : (-rx:rx))
            (dx == 0 && dy == 0) && continue
            nbx, nby = x + dx, y + dy
            (1 <= nbx <= nx && 1 <= nby <= ny) || continue
            j = get(pos_to_idx, nbx + (nby - 1) * nx, 0)
            j == 0 && continue

            ri, rj = find_root(i), find_root(j)
            ri != rj && (parent[ri] = rj)
        end
    end

    # Component sizes come from counting roots - no dense pass.
    roots = Vector{Int}(undef, n)
    comp_size = Dict{Int,Int}()
    for i in 1:n
        r = find_root(i)
        roots[i] = r
        comp_size[r] = get(comp_size, r, 0) + 1
    end

    cleaned_raster = falses(nx, ny)
    for i in 1:n
        comp_size[roots[i]] >= min_cluster_size && (cleaned_raster[coords[i]] = true)
    end

    return cleaned_raster
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
- Masked boolean raster indicating pixels that are within the target distance.
"""
function filter_distances(target_rast::Raster, dist_buffer)::Raster
    # Mask out areas outside considered distance from port
    return mask(Raster(target_rast; missingval=0); with=dist_buffer)
end

"""
    calc_distances(
        target_rast::Raster,
        gdf::DataFrame;
        units::String="NM"
    )::Raster

Calculate the minimum distance from each point in gdf.geometry to each valid pixel in target_rast.

`target_rast` and the `gdf` should be in the same CRS (e.g., EPSG:7844 / GDA2020).

# Arguments
- `target_rast` : Raster of suitable pixels (Bool) to calculate distances from.
- `gdf` : GeoDataFrame of 'points' for distance calculation.

# Returns
- Raster of distances from each cell to the closest point in gdf.
"""
function calc_distances(
    target_rast::Raster,
    gdf::DataFrame;
    units::String="NM"
)::Raster
    tmp_areas = Float32.(copy(target_rast))

    # First dimension is the rows (longitude)
    # Second dimension is the cols (latitude)
    raster_lon = Vector{Float64}(tmp_areas.dims[1].val)
    raster_lat = Vector{Float64}(tmp_areas.dims[2].val)

    # Build KDTree from port coordinates once for O(log n) nearest-port lookup.
    # NearestNeighbors uses Euclidean distance for the tree traversal; Haversine
    # is then computed only for the single nearest port per pixel.
    port_coords = reduce(hcat, [Float64[AG.getx(g, 0), AG.gety(g, 0)] for g in gdf.geometry])
    kdtree = KDTree(port_coords)

    conv = units == "NM" ? 1852.0 : units == "km" ? 1000.0 : 1.0

    valid_idx = findall(!=(0.0f0), tmp_areas.data)

    # Fill a preallocated 2xN matrix directly rather than `reduce(hcat, ...)` over a
    # vector of small per-pixel vectors.
    pixel_coords = Matrix{Float64}(undef, 2, length(valid_idx))
    for (i, idx) in enumerate(valid_idx)
        pixel_coords[1, i] = raster_lon[idx[1]]
        pixel_coords[2, i] = raster_lat[idx[2]]
    end

    # `nn` (not `knn(...,1)`) is NearestNeighbors' dedicated single-nearest-neighbour
    # bulk query - it returns flat `Vector{Int}`/`Vector{Float64}` directly instead of
    # one heap-allocated one-element `Vector{Int}` per query point.
    nearest_idxs, _ = nn(kdtree, pixel_coords)

    # Each iteration writes to a distinct pixel index with no other shared
    # mutable state, so this is safe to run concurrently.
    Threads.@threads for i in eachindex(valid_idx)
        idx = valid_idx[i]
        nearest_port = @view port_coords[:, nearest_idxs[i]]
        dist_nearest = Distances.haversine(nearest_port, @view pixel_coords[:, i])
        tmp_areas.data[idx] = Float32(dist_nearest / conv)
    end

    tmp_areas = rebuild(tmp_areas, missingval=Float32(0.0))
    return tmp_areas
end

"""
    process_UTM_raster(
        src_file::String,
        dst_file::String,
        target_crs::GFT.CoordinateReferenceSystemFormat,
        target_missingval::Float64,
        reg::String,
        method::Symbol
    )::Nothing

Process bathymetry, slope and rugosity datasets from raw input data files and output to `output_fn`
location. These datasets are in region UTM zone CRS.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `src_file` : Path to raw bathymetry or slope raster files in UTM CRS for processing.
- `dst_file` : File location name to create output file. Should include variable and region information.
- `target_crs` : Target CRS object to use in Rasters.resample(). e.g. using GFT.EPSG() format.
- `target_missingval` : Consistent missingval to use in output raster.
- `reg` : Region name for input CRS definition.
- `method` : Resampling interpolation method (more information https://rafaqz.github.io/Rasters.jl/stable/api#Rasters.resample-Tuple).
"""
function process_UTM_raster(
    src_file::String,
    dst_file::String,
    target_crs::GFT.CoordinateReferenceSystemFormat,
    target_missingval::Float64,
    reg::String;
    method::Symbol
)::Nothing
    return skip_if_exists(dst_file) do
        input_raster = Raster(src_file; crs=REGION_CRS_UTM[reg], mappedcrs=EPSG_4326, lazy=true)
        input_raster = set_consistent_missingval!(Rasters.trim(input_raster), target_missingval)

        # GDAL streams block-by-block when writing to disk (what `resample_to_disk`
        # does), so this is the normal path - not a fallback for the rare OOM case.
        resample_to_disk(input_raster, dst_file; crs=target_crs, method=method)

        input_raster = nothing
        force_gc_cleanup()

        return nothing
    end
end

"""
    crop_to_region(
        src_file::String,
        target_region_geom::AbstractVector,
        dst_file::String
    )::Union{Raster,Nothing}

Crop larger input raster to the extent of `target_region_geom` geometry.

# Arguments
- `src_file` : Location of raw input raster file for processing (intended for GBR-wide/rugosity files).
- `target_region_geom` : Region geometry object to crop to (any GeoInterface-compatible geometry vector,
  e.g. `Vector{AG.IGeometry{AG.wkbMultiPolygon}}` or a `GeoDataFrames.GeometryVector`).
- `dst_file` : Path to output file. (File not created within this function, used to check if file already exists).

# Returns
- Raster with the spatial extent matching region_geom.
"""
function crop_to_region(
    src_file::String,
    target_region_geom::AbstractVector,
    dst_file::String
)::Union{Raster,Nothing}
    return skip_if_exists(dst_file) do
        input_raster = Raster(src_file; lazy=true)

        # Note: trim/mask is very important - otherwise file sizes are GBs!
        return Rasters.trim(
            Rasters.mask(
                crop(input_raster; to=target_region_geom); with=target_region_geom
            )
        )
    end
end

"""
    resample_and_write(
        input_raster::Union{Raster,Nothing},
        rst_template::Raster,
        dst_file::String;
        method::Symbol=:near
    )::Nothing

Resample `input_raster` to `template_raster` to ensure matching spatial extent, CRS and resolution.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `input_raster` : Input raster dataset for resampling to template_raster.
- `rst_template` : Template raster for resampling.
- `dst_file` : File location to check - if exists, this function does nothing.
- `method` : Resampling interpolation method supported by Rasters.jl, defaulting to `:near` (nearest neighbor)
             (See [Rasters.jl documentation](https://rafaqz.github.io/Rasters.jl/v0.14.4/api#Rasters.resample-Tuple)).
"""
function resample_and_write(
    input_raster::Union{Raster,Nothing},
    rst_template::Raster,
    dst_file::String;
    method::Symbol=:near
)::Nothing
    return skip_if_exists(dst_file) do
        resample_to_disk(input_raster, dst_file; to=rst_template, method=method)

        return nothing
    end
end

"""
    process_wave_data(
        src_file::String,
        dst_file::String,
        data_layer::Symbol,
        rst_template::Raster,
        target_rst::Raster,
        target_missingval::Float64;
        method::Symbol
    )::Nothing

Process wave data from one CRS/PCS to another, writing the results out to disk as geotiff.

The wave data this function is intended for is provided in netCDF format. The spatial
extents/coordinates provided in this dataset are not well-read by GDAL. To work around this
issue, we use a `rst_template` for the same spatial region to provide replace coordinates,
assuming they are well aligned.

# Notes
- The data is temporarily converted into a sparse matrix to reduce memory use.
  The source data *must not* have valid zero values.
- Existing files are *not* overwritten.

# References
1. Callaghan, David (2023). Great Barrier Reef non-cyclonic and on-reef wave model predictions.
   The University of Queensland.
   Data Collection.
   https://doi.org/10.48610/8246441
   https://espace.library.uq.edu.au/view/UQ:8246441

# Arguments
- `src_file` : Path to netcdf file to process
- `dst_file` : Location of file to write to
- `data_layer` : Name of layer to load
- `rst_template` : Raster in the target Template to use to aid in resampling/reprojection
- `target_rst` : Raster indicating the spatial extent to resample `src_file` into
- `target_missingval` : Intended missingval for the Float64 output
- `method` : Resampling interpolation method (more information https://rafaqz.github.io/Rasters.jl/stable/api#Rasters.resample-Tuple).

# Returns
Nothing
"""
function process_wave_data(
    src_file::String,
    dst_file::String,
    data_layer::Symbol,
    rst_template::Raster,
    target_rst::Raster,
    target_missingval::Float64;
    method::Symbol
)::Nothing
    return skip_if_exists(dst_file; label="Wave data") do
        # Source netCDF files are too large to safely load into memory in one go (the
        # largest region - Mackay-Capricorn - is ~2.9 billion pixels, ~23GB as a plain
        # Julia array), so this stays disk-backed/streamed end to end: crop -> stream to
        # a temp GeoTIFF -> GDAL warp straight to the target grid. Never materialise the
        # full region as a Julia array, and never lazily `reverse` a disk-backed raster
        # (DiskArrays' `readblock!` can't stream a view combining a negative-step
        # `reverse` with crop/broadcast: `MethodError: no method matching readblock!(...,
        # ::StepRange, ...)`) - GDAL's warp handles the source's south-up storage and
        # extending to `target_rst`'s extent natively as part of the resample below.
        #
        # NCDatasets/Rasters doesn't auto-detect this netCDF's `missing_value` attribute
        # as a proper `missingval`, so the raw fill sentinel must be read and passed in
        # explicitly - otherwise it flows through crop/resample as real (wildly
        # out-of-range) data instead of being treated as no-data.
        src_missingval = NCDatasets.NCDataset(src_file) do ds
            Float32(ds[String(data_layer)].attrib["missing_value"])
        end

        wave_rst = Raster(
            src_file,
            name=data_layer,
            crs=GI.crs(rst_template),
            missingval=src_missingval,
            lazy=true
        )

        # Crop to the region's extent - lazy, a plain windowed disk read.
        wave_rst = crop(wave_rst; to=rst_template)
        wave_rst = set_consistent_missingval!(wave_rst, target_missingval)

        tmp_file = tempname() * ".tif"
        try
            Rasters.write(tmp_file, wave_rst; force=true, missingval=Float32(target_missingval))
            wave_rst = nothing
            force_gc_cleanup()

            # Reproject straight to the final GDA2020 grid. GDAL streams block-by-block
            # when writing to disk (what `resample_to_disk` does).
            target_waves = Raster(tmp_file; lazy=true)
            resample_to_disk(target_waves, dst_file; to=target_rst, method=method)
        finally
            isfile(tmp_file) && rm(tmp_file)
        end
        force_gc_cleanup()

        return nothing
    end
end

"""
    distance_raster(
        src_file::String,
        distance_buffer::DataFrame,
        distance_points::DataFrame,
        target_missingval,
        dst_file::String,
        units::String
    )::Nothing

Process and create a raster file containing distance values for each cell from a target geometry.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `src_file` : Path to raster file for processing. Distance will be calculated for all valid pixels.
- `distance_buffer` : DataFrame containing buffer polygon geometries for masking to remove pixels outside of target distance.
- `distance_points` : DataFrame containing point geometries for calculating distance to each pixel. Distance is returned for closest object in distance_points.
- `target_missingval` : Defined missing value to add to output raster for consistency (-9999.0 for GBR-reef-guidance-assessment).
- `dst_file` : Path to output distance .tif file. Should include variable and region information.
- `units` : String of units for output ("m", "km" and "NM" are currently valid inputs).
"""
function distance_raster(
    src_file::String,
    distance_buffer::DataFrame,
    distance_points::DataFrame,
    target_missingval,
    dst_file::String,
    units::String
)::Nothing
    return skip_if_exists(dst_file) do
        target_raster = Raster(src_file; crs=EPSG_7844)
        target_raster = filter_distances(target_raster, distance_buffer)
        target_raster = calc_distances(target_raster, distance_points; units=units)

        target_raster = set_consistent_missingval!(target_raster, target_missingval)
        Rasters.write(dst_file, target_raster)
        target_raster = nothing
        force_gc_cleanup()

        return nothing
    end
end

"""
    within_port_range(
        src_file::String,
        distance_buffer::DataFrame,
        target_missingval::Union{Int64,Float64},
        dst_file::String
    )::Nothing

Create a boolean mask of locations within the distance/range indicated by `distance_buffer`.

Effectively the same method as `distance_raster`, but skipping the distance calculation
step.
"""
function within_port_range(
    src_file::String,
    dist_buffer::DataFrame,
    dst_file::String
)::Nothing
    return skip_if_exists(dst_file) do
        target_raster = Raster(src_file; crs=EPSG_7844, lazy=true)
        target_raster = boolmask(mask(target_raster; with=dist_buffer); missingval=0)
        Rasters.write(dst_file, UInt8.(target_raster))

        target_raster = nothing
        force_gc_cleanup()

        return nothing
    end
end

"""
    write_valid_locs(
        criteria_paths::NamedTuple,
        benthic_ids::Vector,
        geomorph_ids::Vector,
        first_min_size::Int64,
        first_window::Tuple{Int64,Int64},
        second_min_size::Int64,
        second_window::Tuple{Int64,Int64},
        dst_file::String
    )::Nothing

Find the pixels that are covered by valid data for all criteria and benthic/geomorphic IDs.
Applies a two-pass process to remove orphaned pixels - single points of data unconnected to
any other pixel.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `criteria_paths` : NamedTuple containing criteria_fn as keys and path_to_raster_input as values.
- `benthic_ids` : Vector containing the IDs used to filter desired benthic criteria.
- `geomorph_ids` : Vector containing the IDs used to filder desired geomorphic critera (either flat or slope IDs).
- `first_min_size` : Size of minimum cluster to use in `remove_orphaned_elements()` raster cleaning.
- `first_window` : Tuple containing the size of the first window used in cleaning orphaned elements.
- `second_min_size` : Size of minimum cluster to use in `remove_orphaned_elements()` raster cleaning.
- `second_window` : Tuple containing the size of the second window used in cleaning orphaned elements.
- `dst_file` : Path to output results `.tif` file.
"""
function write_valid_locs(
    criteria_paths::NamedTuple,
    benthic_ids::Vector,
    geomorph_ids::Vector,
    first_min_size::Int64,
    first_window::Tuple{Int64,Int64},
    second_min_size::Int64,
    second_window::Tuple{Int64,Int64},
    dst_file::String;
    band_height::Int=2048
)::Nothing
    @assert band_height % 256 == 0 "band_height must be a multiple of 256 to stay aligned with COG blocks"

    return skip_if_exists(dst_file) do
        crits = keys(criteria_paths)
        rsts = NamedTuple{crits}(Tuple(Raster(criteria_paths[c]; lazy=true) for c in crits))
        nx, ny = size(rsts[first(crits)])

        # AND-reduce across all criteria layers, one row band at a time, so no
        # layer's full multi-GB extent is ever materialised - only the current band
        # of each of the 12 layers is resident at once.
        valid_areas = falses(nx, ny)
        y_start = 1
        while y_start <= ny
            y_end = min(y_start + band_height - 1, ny)
            band_acc = nothing
            for crit in crits
                # Slicing the lazy raster keeps this a windowed GDAL read; `boolmask`
                # is applied identically to every criterion (including the class
                # layers, after their `.∈` filter) to exactly replicate the original
                # per-layer missing-value semantics rather than reimplementing them.
                rast = rsts[crit][:, y_start:y_end]
                if crit == :Benthic
                    rast = rast .∈ [benthic_ids]
                elseif crit == :Geomorphic
                    rast = rast .∈ [geomorph_ids]
                end
                band_bool = boolmask(rast)[:, :].data

                if isnothing(band_acc)
                    band_acc = BitMatrix(band_bool)
                else
                    band_acc .&= band_bool
                end
            end
            valid_areas[:, y_start:y_end] .= band_acc
            y_start = y_end + 1
            force_gc_cleanup()
        end

        # Clean up orphaned pixels (first and second pass)
        cleaned_areas = remove_orphaned_elements(BitMatrix(valid_areas), first_min_size, first_window)
        cleaned_areas = remove_orphaned_elements(cleaned_areas, second_min_size, second_window)
        valid_areas = nothing

        # Two-pass trim: scan the BitMatrix for the min/max valid row/col, then
        # build and write UInt8 only over that cropped extent, rather than
        # materialising the full-region UInt8 output (up to 1.24 GB) as
        # `convert.(UInt8, Rasters.trim(...))` on the full raster would.
        x_any = vec(any(cleaned_areas; dims=2))
        y_any = vec(any(cleaned_areas; dims=1))
        x_range = findfirst(x_any):findlast(x_any)
        y_range = findfirst(y_any):findlast(y_any)

        template = rsts[first(crits)][x_range, y_range]
        cropped = Raster(template; data=UInt8.(cleaned_areas[x_range, y_range]), missingval=UInt8(0))

        Rasters.write(dst_file, cropped)

        cleaned_areas = nothing
        force_gc_cleanup()

        return nothing
    end
end

"""
    resize_to_valid_area(criteria_paths::NamedTuple, valid_fn::String)

Resize processed data files to the area that has data across all criteria layers.
Replaces existing file.

# Arguments
- `criteria_paths` : Named collection of criteria paths
- `valid_fn` : Path/name to valid dataset
"""
function resize_to_valid_area(criteria_paths::NamedTuple, valid_fn::String)
    crits = keys(criteria_paths)
    valid_size = size(Raster(valid_fn; lazy=true))

    # Replace files with copies that only cover the relevant valid areas.
    # Each iteration opens its own `valid_areas` handle rather than sharing one
    # opened before the loop, satisfying GDAL's "distinct instance per thread"
    # rule (RFC 101) and making this safe to run under `Threads.@threads`.
    Threads.@threads for crit in collect(crits)
        crit_area = Raster(criteria_paths[crit]; lazy=true)

        if crit ∈ [:Benthic, :Geomorphic]
            m = :near
        else
            m = :bilinear
        end

        if size(crit_area) == valid_size
            @debug "Skipping resizing of $crit ..."
            continue
        end

        valid_areas = Raster(valid_fn; lazy=true)
        trimmed_crit_area = Rasters.trim(crit_area)

        # Disk-based resample (low memory, GDAL streams block-by-block) is the
        # normal path here, not a fallback - the in-memory route this replaced was
        # the only reason these ever OOM'd. This also changes the on-disk encoding
        # of resized criteria layers from striped-plain to 256x256 COG, which is
        # intentional: it's what keeps the windowed reads in `stack_values` and
        # `write_valid_locs` cheap.
        resample_to_disk(trimmed_crit_area, criteria_paths[crit]; to=valid_areas, method=m)
    end

    return nothing
end
