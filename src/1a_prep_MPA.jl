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
    process_wave_data(src_file, dst_file, rst_template, target_rst)::Nothing

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

# Returns
Nothing
"""
function process_wave_data(
    src_file::String,
    dst_file::String,
    data_layer::Symbol,
    rst_template::Raster,
    target_rst::Raster
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
    wave_rst.data[wave_rst.data .< -9999.0] .= -9999.0
    wave_rst = Raster(
        wave_rst;
        data=Float32.(wave_rst.data[:, end:-1:1]),
        missingval=-9999.0
    )

    wave_rst = crop(wave_rst; to=rst_template)

    # Extend bounds of wave data to match bathymetry if needed
    # This is needed to ensure a smaller raster matches the size of the larger raster.
    if !all(size(rst_template) .== size(wave_rst))
        wave_rst = extend_to(wave_rst, rst_template)
        @assert all(size(rst_template) .== size(wave_rst))
    end

    target_Hs = Raster(
        rst_template;
        data=wave_rst.data,
        missingval=-9999.0
    )
    wave_rst = nothing
    force_gc_cleanup()

    # Reproject raster to GDA2020 (degree projection)
    resample(target_Hs; to=target_rst, filename=dst_file, format="COG")
    force_gc_cleanup(; wait_time=2)

    return nothing
end

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    base_bathy_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif")
    if !isfile(base_bathy_fn)
        target_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        target_bathy = Raster(target_bathy_path; mappedcrs=EPSG_4326)
        target_bathy = set_consistent_missingval!(target_bathy, -9999.0)

        resample(target_bathy; crs=EPSG_7844, format="COG", filename=base_bathy_fn)

        target_bathy = nothing
        force_gc_cleanup()
    end

    # Load bathymetry data to provide corresponding spatial extent
    bathy_gda2020 = Raster(base_bathy_fn; lazy=true)

    slope_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif")
    if !isfile(slope_fn)
        target_slope_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        target_slope = Raster(target_slope_path, mappedcrs=EPSG_4326)
        target_slope = set_consistent_missingval!(target_slope, -9999.0)

        resample(target_slope; crs=EPSG_7844, format="COG", filename=slope_fn)

        target_slope = nothing
        force_gc_cleanup(; wait_time=2)
    end

    benthic_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif")
    if !isfile(benthic_fn)
        gbr_benthic_path = "$(MPA_DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
        gbr_benthic = Raster(gbr_benthic_path, mappedcrs=EPSG_4326)

        target_benthic = Rasters.trim(
            mask(
                crop(gbr_benthic; to=regions_4326[reg_idx_4326, :]);
                with=regions_4326[reg_idx_4326, :]
            )
        )

        resample(target_benthic; to=bathy_gda2020, format="COG", filename=benthic_fn)

        gbr_benthic = nothing
        target_benthic = nothing
        force_gc_cleanup(; wait_time=2)
    end

    geomorphic_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif")
    if !isfile(geomorphic_fn)
        gbr_morphic_path = "$(MPA_DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
        target_geomorphic = Raster(gbr_morphic_path, mappedcrs=EPSG_4326)

        target_geomorphic = Rasters.trim(
            mask(
                crop(target_geomorphic; to=regions_4326[reg_idx_4326, :]);
                with=regions_4326[reg_idx_4326, :]
            )
        )

        resample(target_geomorphic; to=bathy_gda2020, format="COG", filename=geomorphic_fn)

        target_geomorphic = nothing
        force_gc_cleanup(; wait_time=2)
    end

    # Use bathy dataset as a template for writing netCDF data to geotiff
    src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    rst_template = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

    waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))
    target_waves_Hs_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif")
    process_wave_data(waves_Hs_path, target_waves_Hs_fn, :Hs90, rst_template, bathy_gda2020)

    waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Tp", reg)))
    target_waves_Tp_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif")
    process_wave_data(waves_Tp_path, target_waves_Tp_fn, :Tp90, rst_template, bathy_gda2020)

    turbid_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif")
    if !isfile(turbid_fn)
        aca_turbid_path = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
        target_turbid = Raster(aca_turbid_path, mappedcrs=EPSG_4326, lazy=true)

        target_turbid = Rasters.trim(
            mask(
                crop(target_turbid; to=bathy_gda2020);
                with=regions_4326[reg_idx_4326, :]
            )
        )

        resample(target_turbid; to=bathy_gda2020, format="COG", filename=turbid_fn)

        target_turbid = nothing
        force_gc_cleanup(; wait_time=2)
    end

    if reg == "Townsville-Whitsunday"
        rugosity_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif")
        if !isfile(rugosity_fn)
            tsv_rugosity_path = joinpath(RUG_DATA_DIR, "std25_Rugosity_Townsville-Whitsunday.tif")
            tsv_rugosity = Raster(tsv_rugosity_path, mappedcrs=EPSG_4326)
            tsv_rugosity = set_consistent_missingval!(tsv_rugosity, -9999.0)

            resample(tsv_rugosity; to=bathy_gda2020, format="COG", filename=rugosity_fn)

            tsv_rugosity = nothing
            force_gc_cleanup(; wait_time=2)
        end
    end

    valid_areas_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif")
    if !isfile(valid_areas_fn)
        # Build validity map/mask
        # Indicates where there is valid data across all criteria

        # Load required prepared raster files for analysis
        rst_template = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"))
        bathy_crit = boolmask(rst_template)
        rst_template = nothing
        force_gc_cleanup(; wait_time=10)  # Needs extra time to clear it seems

        src_slope = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"))
        slope_crit = boolmask(src_slope)
        src_slope = nothing
        force_gc_cleanup(; wait_time=2)

        src_benthic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"))
        benthic_crit = src_benthic .∈ [MPA_BENTHIC_IDS]
        src_benthic = nothing
        force_gc_cleanup(; wait_time=2)

        src_geomorphic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"))
        geomorphic_flat_crit = src_geomorphic .∈ [MPA_FLAT_IDS]
        geomorphic_slope_crit = src_geomorphic .∈ [MPA_SLOPE_IDS]
        src_geomorphic = nothing
        force_gc_cleanup(; wait_time=2)

        src_waves_Hs = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"))
        Hs_waves_crit = boolmask(src_waves_Hs)
        src_waves_Hs = nothing
        force_gc_cleanup(; wait_time=2)

        src_waves_Tp = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"))
        Tp_waves_crit = boolmask(src_waves_Tp)
        src_waves_Tp = nothing
        force_gc_cleanup(; wait_time=2)

        src_turbid = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"))
        turbid_crit = boolmask(src_turbid)
        src_turbid = nothing
        force_gc_cleanup(; wait_time=2)

        # Build mask indicating locations with valid data across all criteria
        valid_areas = (
            benthic_crit .&
            bathy_crit .&
            slope_crit .&
            Hs_waves_crit .&
            Tp_waves_crit .&
            turbid_crit
        )

        if reg == "Townsville-Whitsunday"
            src_rugosity = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"))
            valid_areas .= valid_areas .& boolmask(src_rugosity)
            src_rugosity = nothing
            force_gc_cleanup(; wait_time=2)
        end

        valid_slopes = valid_areas .& geomorphic_slope_crit
        valid_flats = valid_areas .& geomorphic_flat_crit

        # Clean up orphaned pixels (first and second pass)
        cleaned_slopes = remove_orphaned_elements(BitMatrix(valid_slopes.data), 7, (3,3))
        cleaned_slopes = remove_orphaned_elements(cleaned_slopes, 70, (9, 9))
        valid_slopes.data .= cleaned_slopes

        cleaned_flats = remove_orphaned_elements(BitMatrix(valid_flats.data), 7, (3,3))
        cleaned_flats = remove_orphaned_elements(cleaned_flats, 70, (9, 9))
        valid_flats.data .= cleaned_flats

        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif"), convert.(UInt8, valid_slopes))
        write(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif"), convert.(UInt8, valid_flats))

        valid_areas = nothing
        valid_slopes = nothing
        valid_flats = nothing
        cleaned_flats = nothing
        cleaned_slopes = nothing
        force_gc_cleanup(; wait_time=2)
    end

    # Calculate distance to nearest port
    port_dist_slopes_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_slopes.tif")
    if !isfile(port_dist_slopes_fn)
        local port_buffer = GDF.read(joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg"))
        port_points = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))

        valid_slopes = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif"))
        valid_slopes = filter_distances(valid_slopes, port_buffer)
        slope_distances = calc_distances(valid_slopes, port_points; units="NM")

        slope_distances = set_consistent_missingval!(slope_distances, -9999.0)
        write(port_dist_slopes_fn, slope_distances)
        valid_slopes = nothing
        slope_distances = nothing
        force_gc_cleanup(; wait_time=2)
    end

    port_dist_flats_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_flats.tif")
    if !isfile(port_dist_flats_fn)
        local port_buffer = GDF.read(joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg"))
        port_points = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))

        valid_flats = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif"))
        valid_flats = filter_distances(valid_flats, port_buffer)
        flat_distances = calc_distances(valid_flats, port_points; units="NM")

        flat_distances = set_consistent_missingval!(flat_distances, -9999.0)
        write(port_dist_flats_fn, flat_distances)
        valid_flats = nothing
        flat_distances = nothing
        force_gc_cleanup(; wait_time=2)
    end

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
