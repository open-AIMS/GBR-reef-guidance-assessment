using
    Rasters,
    ImageMorphology,
    DataFrames

import ArchGDAL as AG
import GeoInterface as GI
import GeoFormatTypes as GFT


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
    process_bottom_UTM_raster(
        src_file::String,
        dst_file::String,
        target_crs::GFT.CoordinateReferenceSystemFormat,
        target_missingval::Float64,
        reg::String
    )::Nothing

Process bathymetry, slope and rugosity datasets from raw input data files and output to `output_fn`
location. These datasets are in region UTM zone CRS.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `src_file` : Path to raw bathymetry, slope or rugosity raster files in UTM CRS for processing.
- `dst_file` : File location name to create output file. Should include variable and region information.
- `target_crs` : Target CRS object to use in Rasters.resample(). e.g. using GFT.EPSG() format.
- `target_missingval` : Consistent missingval to use in output raster.
- `reg` : Region name for input CRS definition.

"""
function process_bottom_UTM_raster(
    src_file::String,
    dst_file::String,
    target_crs::GFT.CoordinateReferenceSystemFormat,
    target_missingval::Float64,
    reg::String
)::Nothing
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    input_raster = Raster(src_file; crs=REGION_CRS_UTM[reg], mappedcrs=EPSG_4326)
    input_raster = set_consistent_missingval!(input_raster, target_missingval)

    resample(input_raster; crs=target_crs, filename=dst_file, format="COG")
    force_gc_cleanup(; wait_time=2)

    return nothing
end

"""
    trim_extent_region(
        src_file::String,
        input_crs::GFT.CoordinateReferenceSystemFormat,
        target_region_geom::Vector{AG.IGeometry{AG.wkbMultiPolygon}},
        dst_file::String
    )::Union{Raster,Nothing}

Trim larger input raster to the extent of `region_geom` geometry.

# Arguments
- `src_file` : Location of raw input raster file for processing (intended for GBR-wide/rugosity files).
- `input_crs` : CRS of input raster file using GFT.EPSG().
- `target_region_geom` : Region geometry object to crop to.
- `dst_file` : Path to output file. (File not created within this function, used to check if file already exists).

# Returns
- Raster with the spatial extent matching region_geom.
"""
function trim_extent_region(
    src_file::String,
    input_crs::GFT.CoordinateReferenceSystemFormat,
    target_region_geom::Vector{AG.IGeometry{AG.wkbMultiPolygon}},
    dst_file::String
)::Union{Raster,Nothing}
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    input_raster = Raster(src_file; mappedcrs=input_crs, lazy=true)

    return Rasters.trim(
        mask(
            crop(input_raster; to=target_region_geom);
            with=target_region_geom
        )
    )
end

"""
    resample_and_write(
        input_raster::Union{Raster,Nothing},
        template_raster::Raster
        dst_file::String
    )::Nothing

Resample `input_raster` to `template_raster` to ensure matching spatial extent, CRS and resolution.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `input_raster` : Input raster dataset for resampling to template_raster.
- `rst_template` : Template raster for resampling.
- `dst_file` : File location name to create output file. Should include variable and region information.
"""
function resample_and_write(
    input_raster::Union{Raster,Nothing},
    rst_template::Raster,
    dst_file::String
)::Nothing
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    resample(input_raster; to=rst_template, filename=dst_file, format="COG")

    return nothing
end

"""
    process_wave_data(
        src_file::String,
        dst_file::String,
        data_layer::Symbol,
        rst_template::Raster,
        target_rst::Raster,
        target_missingval::Float64
    )::Nothing

Process wave data from one CRS/PCS to another, writing the results out to disk in COG
format.

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

# Returns
Nothing
"""
function process_wave_data(
    src_file::String,
    dst_file::String,
    data_layer::Symbol,
    rst_template::Raster,
    target_rst::Raster,
    target_missingval::Float64
)::Nothing
    if isfile(dst_file)
        @warn "Wave data not processed as $(dst_file) already exists."
        return
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
    wave_rst.data[wave_rst.data .< target_missingval] .= target_missingval
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
    resample(target_waves; to=target_rst, filename=dst_file, format="COG")
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
- `target_missingval` : missingval to add to output raster for consistency (-9999.0 for GBR-reef-guidance-assessment).
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
    write(dst_file, target_raster)
    target_raster = nothing
    force_gc_cleanup()

    return nothing
end

"""
    find_valid_locs(
        criteria_paths::NamedTuple,
        benthic_ids::Vector,
        geomorph_ids::Vector,
        first_min_size::Int64,
        first_window::Tuple{Int64, Int64},
        second_min_size::Int64,
        second_window::Tuple{Int64, Int64},
        dst_file::String,
        reg::String
    )::Nothing

Find the pixels that are covered by valid data for all criteria and benthic/geomorphic IDs.

Writes to `dst_file` as a Cloud Optimized Geotiff.

# Arguments
- `criteria_paths` : NamedTuple containing criteria_fn as keys and path_to_raster_input as values.
- `benthic_ids` : Vector containing the IDs used to filter desired benthic criteria.
- `geomorph_ids` : Vector containing the IDs used to filder desired geomorphic critera (either flat or slope IDs).
- `first_min_size` : Size of minimum cluster to use in `remove_orphaned_elements()` raster cleaning.
- `first_window` : Tuple containing the size of the first window used in cleaning orphaned elements.
- `second_min_size` : Size of minimum cluster to use in `remove_orphaned_elements()` raster cleaning.
- `second_window` : Tuple containing the size of the second window used in cleaning orphaned elements.
- `dst_file` : Path to output results .tif file.
- `reg` : Current processing region.
"""
function find_valid_locs(
    criteria_paths::NamedTuple,
    benthic_ids::Vector,
    geomorph_ids::Vector,
    first_min_size::Int64,
    first_window::Tuple{Int64,Int64},
    second_min_size::Int64,
    second_window::Tuple{Int64,Int64},
    dst_file::String,
    reg::String
)::Nothing
    if isfile(dst_file)
        @warn "Data not processed as $(dst_file) already exists."
        return
    end

    src_bathy = Raster(criteria_paths[:bathy_fn])
    bathy_crit = boolmask(src_bathy)
    rst_template = nothing
    force_gc_cleanup(; wait_time=10)  # Needs extra time to clear it seems

    src_slope = Raster(criteria_paths[:slope_fn])
    slope_crit = boolmask(src_slope)
    src_slope = nothing
    force_gc_cleanup(; wait_time=2)

    src_benthic = Raster(criteria_paths[:benthic_fn])
    benthic_crit = src_benthic .∈ [benthic_ids]
    src_benthic = nothing
    force_gc_cleanup(; wait_time=2)

    src_geomorphic = Raster(criteria_paths[:geomorph_fn])
    geomorphic_crit = src_geomorphic .∈ [geomorph_ids]
    src_geomorphic = nothing
    force_gc_cleanup(; wait_time=2)

    src_waves_Hs = Raster(criteria_paths[:waves_Hs_fn])
    Hs_waves_crit = boolmask(src_waves_Hs)
    src_waves_Hs = nothing
    force_gc_cleanup(; wait_time=2)

    src_waves_Tp = Raster(criteria_paths[:waves_Tp_fn])
    Tp_waves_crit = boolmask(src_waves_Tp)
    src_waves_Tp = nothing
    force_gc_cleanup(; wait_time=2)

    src_turbid = Raster(criteria_paths[:turbid_fn])
    turbid_crit = boolmask(src_turbid)
    src_turbid = nothing
    force_gc_cleanup(; wait_time=2)

    # Build mask indicating locations with valid data across all criteria
    valid_areas = (
        benthic_crit .&
        bathy_crit .&
        slope_crit .&
        Hs_waves_crit .&
        turbid_crit .&
        Tp_waves_crit .&
        geomorphic_crit
    )

    if reg == "Townsville-Whitsunday"
        src_rugosity = Raster(criteria_paths[:rugosity_fn])
        valid_areas .= valid_areas .& boolmask(src_rugosity)
        src_rugosity = nothing
        force_gc_cleanup(; wait_time=2)
    end

    # Clean up orphaned pixels (first and second pass)
    cleaned_areas = remove_orphaned_elements(BitMatrix(valid_areas.data), first_min_size, first_window)
    cleaned_areas = remove_orphaned_elements(cleaned_areas, second_min_size, second_window)
    valid_areas.data .= cleaned_areas

    write(dst_file, convert.(UInt8, valid_areas))

    valid_areas = nothing
    cleaned_areas = nothing
    force_gc_cleanup(; wait_time=2)

    return nothing
end
