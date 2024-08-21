"""
Draft raster and polygon processing functions for ease of use in data-preprocessing.
"""

"""
        function process_MPA_bottom_raster(
        input_raster,
        target_crs,
        target_missingval,
        output_fn
    )

Process MPA bathymetry and slope datasets from raw input data files and output to `output_fn`
location.

# Arguments
- `input_raster` : Raw GBRMPA bathymetry or slope raster file
- `target_crs` : Target CRS object to use in Rasters.resample(). e.g. using GFT.EPSG() format.
- `target_missingval` : Consistent missingval to use in output raster.
- `output_fn` : File location name to create output file. Should include variable and region information.
"""
function process_MPA_bottom_raster(
    input_raster::Raster,
    target_crs,
    target_missingval,
    output_fn::String
)
    if !isfile(output_fn)
        input_raster = set_consistent_missingval(input_raster, target_missingval)

        write(
            output_fn,
            resample(
                input_raster;
                crs=target_crs
            )
        )

        input_raster = nothing
        force_gc_cleanup()
    end
end

"""
        function process_MPA_bottom_raster(
        input_raster_fn,
        input_crs,
        target_crs,
        target_missingval,
        output_fn
    )

Process MPA bathymetry and slope datasets from raw input data files and output to `output_fn`
location.

# Arguments
- `input_raster_fn` : Location of raw GBRMPA bathymetry or slope raster file for processing.
- `input_crs` : CRS of input raster file in GFT.EPSG format.
- `target_crs` : Target CRS object to use in Rasters.resample(). e.g. using GFT.EPSG() format.
- `target_missingval` : Consistent missingval to use in output raster.
- `output_fn` : File location name to create output file. Should include variable and region information.
"""
function process_MPA_bottom_raster(
    input_raster_fn,
    input_crs,
    target_crs,
    target_missingval,
    output_fn
)
    if !isfile(output_fn)
        input_raster = Raster(input_raster_fn; mappedcrs=input_crs)
        input_raster = set_consistent_missingval(input_raster, target_missingval)

        write(
            output_fn,
            resample(
                input_raster;
                crs=target_crs
            )
        )

        input_raster = nothing
        force_gc_cleanup()
    end
end

"""
        function trim_extent_region(
        input_raster_fn,
        input_crs,
        region_geom
    )

Trim larger input raster to the extent of region_geom geometry.

# Arguments
- `input_raster_fn` : Location of raw input raster file for processing (intended for GBR-wide/rugosity files).
- `input_crs` : CRS of input raster file in GFT.EPSG format.
- `region_geom` : Region geometry object to crop to.

# Returns
- Raster with the spatial extent matching region_geom.
"""
function trim_extent_region(
    input_raster_fn,
    input_crs,
    region_geom
)
    input_raster = Raster(input_raster_fn; mappedcrs=input_crs)

    return Rasters.trim(
        mask(
            crop(input_raster; to=region_geom);
            with=region_geom
        )
    )
end

"""
        function resample_and_write(
        input_raster,
        template_raster
        output_fn
    )

Resample input_raster to template_raster to ensure matching spatial extent, CRS and resolution.

# Arguments
- `input_raster` : Raster file for resampling to template_raster.
- `template_raster` : Template raster for resampling.
- `output_fn` : File location name to create output file. Should include variable and region information.
"""
function resample_and_write(
    input_raster,
    template_raster,
    output_fn
)
    if !isfile(output_fn)
        write(
            output_fn,
            resample(
                input_raster;
                to=template_raster
            )
        )
    end
end

"""
        function distance_raster(
        template_raster_fn,
        distance_buffer::DataFrame,
        distance_points::DataFrame,
        missingval,
        output_fn,
        units::String
    )

Process and create a raster file containing distance values for each cell from a target geometry.

# Arguments
- `template_bool_fn` : Location of raster file for processing. In `filter_distances()` step the raster will be converted to bool values. Distance will be calculated for all valid pixels.
- `distance_buffer` : DataFrame containing buffer polygon geometries for masking to remove pixels outside of target distance.
- `distance_points` : DataFrame containing point geometries for calculating distance to each pixel. Distance is returned for closest object in distance_points.
- `target_missingval` : missingval to add to output raster for consistency (-9999.0 for GBR-reef-guidance-assessment).
- `output_fn` : Location to output distance .tif file. Should include variable and region information.
- `units` : String of units for output ("m", "km" and "NM" are currently valid inputs).
"""
function distance_raster(
    template_bool_fn,
    distance_buffer::DataFrame,
    distance_points::DataFrame,
    target_missingval,
    output_fn,
    units::String
)
    if !isfile(output_fn)
        target_raster = Raster(template_bool_fn; crs=EPSG_7844)
        target_raster = filter_distances(target_raster, distance_buffer)
        target_raster = calc_distances(target_raster, distance_points; units=units)

        target_raster = set_consistent_missingval(target_raster, missingval)
        write(output_fn, target_raster)
        target_raster = nothing
        force_gc_cleanup()
    end
end
