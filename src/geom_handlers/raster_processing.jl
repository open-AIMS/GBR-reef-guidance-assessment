using
    Rasters,
    ImageMorphology,
    DataFrames

import ArchGDAL as AG
import GeoInterface as GI
import GeoFormatTypes as GFT


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
    replace_missing!(raster, val)
    return Raster(raster; missingval=val)
end

"""
    extend_to(rst1::Raster, rst2::Raster)::Raster

Extend bounds of a `rst1` to the same shape as `rst2`
"""
function extend_to(rst1::Raster, rst2::Raster)::Raster
    rst1 = extend(rst1; to=GI.extent(rst2))
    @assert all(size(rst1) .== size(rst2)) "Sizes do not match post-extension: $(size(rst1)) $(size(rst2))"

    return rst1
end

"""
    remove_orphaned_elements(rst_mask::BitMatrix, min_cluster_size::Int, box_size::Tuple{Int64,Int64})

Cleans up valid pixels that are by themselves and not worth including in later assessments.

# Arguments
- `rst_mask` : Mask of valid raster locations
- `min_cluster_size` : Number of elements that need to be clustered together to be kept
- `box_size` : area to search around center pixel (width, height). Must be odd numbers.
"""
function remove_orphaned_elements(rst_mask::BitMatrix, min_cluster_size::Int, box_size::Tuple{Int64,Int64})
    labels = label_components(rst_mask, strel_box(box_size))

    # Count the size of each component
    component_sizes = component_lengths(labels)

    # Mask components to keep
    keep_mask = component_sizes .>= min_cluster_size

    cleaned_raster = map(x -> keep_mask[x], labels) .* rst_mask

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

    FLoops.assistant(false)
    @floop for (lon_ind, lon) in enumerate(raster_lon)
        for (lat_ind, lat) in enumerate(raster_lat)
            if tmp_areas[lon_ind, lat_ind] != 0.0
                point = AG.createpoint()
                AG.addpoint!(point, lon, lat)

                pixel_dists = AG.distance.([point], gdf.geometry)
                geom_point = gdf[argmin(pixel_dists), :geometry]
                geom_point = (AG.getx(geom_point, 0), AG.gety(geom_point, 0))

                dist_nearest = Distances.haversine(geom_point, (lon, lat))

                # Convert from meters to nautical miles
                if units == "NM"
                    dist_nearest = dist_nearest / 1852
                end

                # Convert from meters to kilometers
                if units == "km"
                    dist_nearest = dist_nearest / 1000
                end

                tmp_areas.data[lon_ind, lat_ind] = Float32(dist_nearest)
            end
        end
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
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return nothing
    end

    input_raster = Raster(src_file; crs=REGION_CRS_UTM[reg], mappedcrs=EPSG_4326)
    input_raster = set_consistent_missingval!(Rasters.trim(input_raster), target_missingval)

    try
        # Using `filename` argument reduces memory use but the resulting file is
        # orders of magnitude bigger, so write out manually
        input_raster = resample(input_raster; crs=target_crs, method=method)
        Rasters.write(dst_file, input_raster)
    catch
        # In cases where data is too big to fit into memory, attempt to resample
        # writing out to disk. This results in a larger than usual file but better
        # than crashing out.
        resample(input_raster; crs=target_crs, method=method, filename=dst_file)
    end

    input_raster = nothing
    force_gc_cleanup(; wait_time=2)

    return nothing
end

"""
    crop_to_region(
        src_file::String,
        target_region_geom::Vector{AG.IGeometry{AG.wkbMultiPolygon}},
        dst_file::String
    )::Union{Raster,Nothing}

Crop larger input raster to the extent of `target_region_geom` geometry.

# Arguments
- `src_file` : Location of raw input raster file for processing (intended for GBR-wide/rugosity files).
- `target_region_geom` : Region geometry object to crop to.
- `dst_file` : Path to output file. (File not created within this function, used to check if file already exists).

# Returns
- Raster with the spatial extent matching region_geom.
"""
function crop_to_region(
    src_file::String,
    target_region_geom::Vector{AG.IGeometry{AG.wkbMultiPolygon}},
    dst_file::String
)::Union{Raster,Nothing}
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return nothing
    end

    input_raster = Raster(src_file; lazy=true)

    # Note: trim/mask is very important - otherwise file sizes are GBs!
    return Rasters.trim(
        Rasters.mask(
            crop(input_raster; to=target_region_geom); with=target_region_geom
        )
    )

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
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return nothing
    end

    # Using `filename` argument reduces memory use but explodes size of file.
    # https://github.com/rafaqz/Rasters.jl/issues/706
    input_raster = resample(
        Raster(input_raster; data=ExtendableSparse.sparse(coalesce.(input_raster.data, 0)));
        to=rst_template,
        method=method
    )
    Rasters.write(dst_file, input_raster)

    return nothing
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
    if isfile(dst_file)
        @warn "Wave data not processed as $(dst_file) already exists."
        return nothing
    end

    # Have to load netCDF data into memory to allow missing value replacement
    wave_rst = Raster(
        src_file,
        name=data_layer,
        crs=GI.crs(rst_template),
        mappedcrs=EPSG_4326
    )

    # 1. Manually set -infinite missing data value to exact value
    #    This is necessary as the netCDF was provided without a set `no data` value
    # 2. We also want to make the type explicit, from Union{Missing,Float32} -> Float32
    # 3. Important to flip the y-axis as the data was stored in reverse orientation
    #    (south-up), so we flip it back (2nd dimension is the y-axis)
    wave_rst.data[wave_rst.data.<target_missingval] .= target_missingval
    wave_rst = Raster(
        wave_rst;
        data=Float32.(wave_rst.data[:, end:-1:1]),
        missingval=target_missingval
    )

    wave_rst = crop(wave_rst; to=rst_template)

    # Extend bounds of wave data to match bathymetry if needed
    # This is needed to ensure a smaller raster matches the size of the larger raster.
    if !all(size(rst_template) .== size(wave_rst))
        wave_rst = extend_to(wave_rst, rst_template)
        @assert all(size(rst_template) .== size(wave_rst))
    end

    target_waves = Raster(
        rst_template;
        data=wave_rst.data,
        missingval=target_missingval
    )
    wave_rst = nothing
    force_gc_cleanup()

    # Reproject raster to GDA2020 (degree projection)
    # Using `filename` argument reduces memory use but the resulting file is
    # orders of magnitude bigger, so write out manually
    target_waves = resample(target_waves; to=target_rst, method=method)
    Rasters.write(dst_file, target_waves)
    force_gc_cleanup(; wait_time=2)

    return nothing
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
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    target_raster = Raster(src_file; crs=EPSG_7844)
    target_raster = filter_distances(target_raster, distance_buffer)
    target_raster = calc_distances(target_raster, distance_points; units=units)

    target_raster = set_consistent_missingval!(target_raster, target_missingval)
    Rasters.write(dst_file, target_raster)
    target_raster = nothing
    force_gc_cleanup()

    return nothing
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
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return nothing
    end

    target_raster = Raster(src_file; crs=EPSG_7844)
    target_raster = boolmask(mask(target_raster; with=dist_buffer); missingval=0)
    Rasters.write(dst_file, UInt8.(target_raster))

    target_raster = nothing
    force_gc_cleanup()

    return nothing
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
    dst_file::String
)::Nothing
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return nothing
    end

    crits = keys(criteria_paths)
    valid_areas = boolmask(Raster(criteria_paths[crits[1]]))
    for crit in crits[2:end]
        rast = Raster(criteria_paths[crit])
        if crit == :Benthic
            rast = rast .∈ [benthic_ids]
        elseif crit == :Geomorphic
            rast = rast .∈ [geomorph_ids]
        end

        valid_areas = valid_areas .& boolmask(rast)
    end

    # Clean up orphaned pixels (first and second pass)
    cleaned_areas = remove_orphaned_elements(BitMatrix(valid_areas.data), first_min_size, first_window)
    cleaned_areas = remove_orphaned_elements(cleaned_areas, second_min_size, second_window)
    valid_areas.data .= cleaned_areas
    valid_areas = convert.(UInt8, Rasters.trim(valid_areas))

    Rasters.write(dst_file, valid_areas)

    cleaned_areas = nothing
    force_gc_cleanup(; wait_time=2)

    return nothing
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
    valid_areas = Raster(valid_fn; lazy=true)

    # Replace files with copies that only cover the relevant valid areas
    for crit in crits
        crit_area = Raster(criteria_paths[crit]; lazy=true)

        if crit ∈ [:Benthic, :Geomorphic]
            m = :near
        else
            m = :bilinear
        end

        if size(crit_area) == size(valid_areas)
            @debug "Skipping resizing of $crit ..."
            continue
        end

        # Using `filename` argument reduces memory use but the resulting file is
        # orders of magnitude bigger, so write out manually
        crit_area = resample(
            read(Rasters.trim(crit_area));
            to=valid_areas,
            method=m
        )
        Rasters.write(criteria_paths[crit], crit_area; force=true)
    end

    return nothing
end
