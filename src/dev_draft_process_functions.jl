"""
Draft raster and polygon processing functions for ease of use in data-preprocessing.
"""

function process_MPA_bottom_raster(
    input_raster,
    target_crs,
    target_missingval,
    output_fn
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

function distance_raster(
    template_raster_fn,
    distance_buffer::DataFrame,
    distance_points::DataFrame,
    missingval,
    output_fn,
    units::String
)
    if !isfile(output_fn)
        target_raster = Raster(template_raster_fn; crs=EPSG_7844)
        target_raster = filter_distances(target_raster, distance_buffer)
        target_raster = calc_distances(target_raster, distance_points; units=units)

        target_raster = set_consistent_missingval(target_raster, missingval)
        write(output_fn, target_raster)
        target_raster = nothing
        force_gc_cleanup()
    end
end
