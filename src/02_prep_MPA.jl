"""
Prepare data for analysis by processing UQ-GBRMPA files for each GBRMPA management region.

"MPA" is the internal designation. This dataset is commonly known as the UQ-GBRMPA or,
for the bathymetry data specifically, the EoMap dataset.

The general approach is to crop down to an extent of the region, then trim/mask to just the
areas of interest. Processing data in this way minimizes the amount of data being handled.

Still, the largest chunk of time is spent resampling to the same projections and writing
data to file (as it takes time to compress the data). The datasets can be in the GBs without
compression.

The steps are:
1. Crop GBR-wide GBRMPA rasters into management regions.
2. Reproject all data from WGS84 / UTM Zone 54 - 56 into a consistent CRS (GDA2020).
3. Ensure all rasters are the same size/shape for each region of interest with the same \
   values used to indicate no data.
"""

include("common.jl")
include("geom_handlers/geom_ops.jl")
include("geom_handlers/raster_processing.jl")
include("geom_handlers/lookup_processing.jl")

using Dates
using SparseArrays, NamedTupleTools
using ExtendableSparse
using Distributed

# Orphaned-pixel cleanup parameters for `write_valid_locs` (first pass, second pass)
const ORPHAN_CLEANUP_FIRST_MIN_SIZE = 7
const ORPHAN_CLEANUP_FIRST_WINDOW = (3, 3)
const ORPHAN_CLEANUP_SECOND_MIN_SIZE = 70
const ORPHAN_CLEANUP_SECOND_WINDOW = (9, 9)

"""
    load_region_geometries()::DataFrame

Load the region boundary geometries (EPSG:4326) used to crop vector and raster
data down to each GBRMPA management region.
"""
function load_region_geometries()::DataFrame
    return GDF.read(REGION_PATH_4326)
end

"""
    prep_ports_and_zoning()::String

Prepare the inputs that are shared across all regions:
1a. Process the GBRMPA zoning geopackage (used to avoid Pink/Exclusion zones).
1b. Reproject port locations into GDA2020.
1c. Build a distance buffer around ports (`PORT_RANGE_NM` is the max range considered).

Returns the path to the written port buffer geopackage rather than the
in-memory geometry: LibGEOS/ArchGDAL objects are not safely shareable across
worker processes, so per-region processing re-reads this file itself instead
of receiving a live handle.
"""
function prep_ports_and_zoning()::String
    # 1a. Process GBRMPA zoning geopackage (used to avoid Pink/Exclusion zones)
    MPA_zoning_input = "$(GDA2020_DATA_DIR)/Great_Barrier_Reef_Marine_Park_Zoning_20_4418126048110066699.gpkg"
    MPA_preservation_zone_fn = joinpath(MPA_OUTPUT_DIR, "GBRMPA_preservation_zone_exclusion.gpkg")
    geometry_exclusion_process(
        MPA_zoning_input,
        MPA_preservation_zone_fn,
        EPSG_7844,
        EPSG_7844,
        :TYPE,
        MPA_EXCLUSION_ZONES;
        geom_col=:SHAPE
    )

    # 1b. Process port locations
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

    # 1c. Create buffer around ports (PORT_RANGE_NM said to be max range considered)
    port_buffer_fn = joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg")
    if !isfile(joinpath(MPA_OUTPUT_DIR, "ports_buffer.gpkg"))
        port_locs = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))

        port_buffer = port_buffer_mask(port_locs, PORT_RANGE_NM, unit="NM")
        port_buffer = DataFrame(; Name="ports_buffer", geometry=port_buffer)
        GDF.write(
            port_buffer_fn,
            port_buffer;
            crs=EPSG_7844
        )
    end

    return port_buffer_fn
end

"""
    write_region_outputs(reg::String, criteria_paths)::Nothing

Final "output write" stage of per-region MPA processing: compute the region's
valid-data mask, resize all criteria rasters down to that valid area, and
write the Arrow lookup table used for fast querying.
"""
function write_region_outputs(reg::String, criteria_paths)::Nothing
    @debug "$(now()) - Processing $(reg) - Valid area"
    # Find locations containing valid data
    valid_slopes_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif")
    if !isfile(valid_slopes_fn)
        write_valid_locs(
            criteria_paths,
            [values(MPA_BENTHIC_IDS)...],
            [values(MPA_SLOPE_IDS)...],
            ORPHAN_CLEANUP_FIRST_MIN_SIZE, ORPHAN_CLEANUP_FIRST_WINDOW,
            ORPHAN_CLEANUP_SECOND_MIN_SIZE, ORPHAN_CLEANUP_SECOND_WINDOW,
            valid_slopes_fn
        )

        resize_to_valid_area(
            criteria_paths,
            valid_slopes_fn
        )
    end

    # Create lookup tables to support fast querying
    @debug "$(now()) - Processing $(reg) - Lookup table"
    slopes_lookup_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes_lookup.arrow")
    valid_lookup(
        criteria_paths,
        valid_slopes_fn,
        slopes_lookup_fn
    )

    return nothing
end

"""
    process_region(reg::String, port_buffer_fn::String)::Nothing

Process a single GBRMPA management region: crop/reproject/resample all
criteria rasters (bathymetry, slope, benthic, geomorphic, turbidity,
rugosity, waves, tides, port distance), then write the region's final
valid-area outputs and lookup table (see `write_region_outputs`).

Only serializable arguments are accepted (the region name and a path to the
port buffer geopackage) so this function can be called from a remote
`Distributed.jl` worker: any GEOS/GDAL state it needs (`regions_4326`,
`port_buffer`) is rebuilt locally from disk rather than passed in as a live
object, since LibGEOS/ArchGDAL contexts are not safe to share across
processes.
"""
function process_region(reg::String, port_buffer_fn::String)::Nothing
    regions_4326 = load_region_geometries()
    reg_idx_4326 = region_geom_index(regions_4326, reg)
    port_buffer = GDF.read(port_buffer_fn)

    # Create NamedTuple to hold all output file paths.
    criteria_paths = create_criteria_paths(reg)

    @debug "$(now()) - Processing $(reg) - bathy and slope"

    # Process bathymetry and slope UTM raster files
    # Write to "[some_file].tif.tif" temporarily.
    # The extra extension is used so the correct format is auto-selected without
    # needing a separate temporary filename (useful for debugging)
    if !isfile(criteria_paths[:Depth])
        raw_bathy_fn = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        process_UTM_raster(
            raw_bathy_fn,
            criteria_paths[:Depth] * ".tif",
            EPSG_7844,
            DEFAULT_MISSINGVAL,
            reg;
            method=:bilinear
        )

        target_depth = crop_to_region(
            criteria_paths[:Depth] * ".tif",
            regions_4326[reg_idx_4326, :geometry],
            criteria_paths[:Depth]
        )

        if !isnothing(target_depth)
            # Write out cropped dataset if needed
            Rasters.write(criteria_paths[:Depth], target_depth; force=true)
            target_depth = nothing
        end

        # Delete the temporary copy
        rm(criteria_paths[:Depth] * ".tif"; force=true)
        force_gc_cleanup()
    end

    # Load bathymetry data to provide corresponding spatial extent
    bathy_gda2020 = Raster(criteria_paths[:Depth]; crs=EPSG_7844, lazy=true)
    # template_rst = Raster(bathy_gda2020; data=zeros(Float32, size(bathy_gda2020)))

    # Write to "[some_file].tif.tif" temporarily.
    # The extra extension is used so the correct format is auto-selected without
    # needing a separate temporary filename (useful for debugging)
    if !isfile(criteria_paths[:Slope])
        raw_slope_fn = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        process_UTM_raster(
            raw_slope_fn,
            criteria_paths[:Slope] * ".tif",
            EPSG_7844,
            DEFAULT_MISSINGVAL,
            reg;
            method=:bilinear
        )

        target_slope = crop_to_region(
            criteria_paths[:Slope] * ".tif",
            regions_4326[reg_idx_4326, :geometry],
            criteria_paths[:Slope]
        )
        force_gc_cleanup()

        # Resample to align with depth dataset
        resample_and_write(
            target_slope,
            bathy_gda2020,
            criteria_paths[:Slope];
            method=:bilinear
        )

        # Delete the temporary copy
        rm(criteria_paths[:Slope] * ".tif"; force=true)

        target_slope = nothing
        force_gc_cleanup()
    end

    # Process other raster data for region
    @debug "$(now()) - Processing $(reg) - Benthic"
    resample_and_write(
        Raster(create_intermediate_filenames(reg).Benthic; lazy=true),
        bathy_gda2020,
        criteria_paths[:Benthic];
        method=:near
    )
    force_gc_cleanup()

    @debug "$(now()) - Processing $(reg) - Geomorphic"
    resample_and_write(
        Raster(create_intermediate_filenames(reg).Geomorphic; lazy=true),
        bathy_gda2020,
        criteria_paths[:Geomorphic];
        method=:near
    )
    force_gc_cleanup()

    @debug "$(now()) - Processing $(reg) - Turbidity"
    if !isfile(criteria_paths[:Turbidity])
        raw_turbid_fn = "$(ACA_DATA_DIR)/Turbidity-2024/turbidity-annual_0.tif"
        target_turbid = crop_to_region(
            raw_turbid_fn,
            regions_4326[reg_idx_4326, :geometry],
            criteria_paths[:Turbidity]
        )
        force_gc_cleanup()

        resample_and_write(
            target_turbid,
            bathy_gda2020,
            criteria_paths[:Turbidity];
            method=:bilinear
        )
        target_turbid = nothing
        force_gc_cleanup()
    end

    # Process Rugosity data
    if reg == "Townsville-Whitsunday"
        @debug "$(now()) - Processing $(reg) - Rugosity"
        raw_rugosity_fn = joinpath(RUG_DATA_DIR, "std25_Rugosity_Townsville-Whitsunday.tif")
        resample_and_write(
            Raster(raw_rugosity_fn; crs=REGION_CRS_UTM[reg], mappedcrs=EPSG_4326, lazy=true),
            bathy_gda2020,
            criteria_paths[:Rugosity];
            method=:bilinear
        )
    end

    # Process wave raster data
    # Use bathy dataset as a template for writing netCDF data to geotiff
    src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    rst_template = Raster(
        src_bathy_path;
        crs=REGION_CRS_UTM[reg],
        mappedcrs=EPSG(4326),
        lazy=true
    )

    @debug "$(now()) - Processing $(reg) - Waves Hs"
    waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))
    process_wave_data(
        waves_Hs_path,
        criteria_paths[:WavesHs],
        :Hs90,
        rst_template,
        bathy_gda2020,
        DEFAULT_MISSINGVAL;
        method=:bilinear
    )

    @debug "$(now()) - Processing $(reg) - Waves Tp"
    waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Tp", reg)))
    process_wave_data(
        waves_Tp_path,
        criteria_paths[:WavesTp],
        :Tp90,
        rst_template,
        bathy_gda2020,
        DEFAULT_MISSINGVAL;
        method=:bilinear
    )

    @debug "$(now()) - Processing $(reg) - Waves Ubed90"
    waves_Ubed_path = first(glob("*ubed90*.nc", joinpath(WAVE_DATA_DIR, "Ubed", reg)))
    process_wave_data(
        waves_Ubed_path,
        criteria_paths[:WavesUbed],
        :ubed90,
        rst_template,
        bathy_gda2020,
        DEFAULT_MISSINGVAL;
        method=:bilinear
    )

    # Tidal data is already separated into management regions
    # so we only need to reproject into consistent datum
    @debug "$(now()) - Processing $(reg) - High Tide"
    hightide_fn = first(glob("*_hightide_*_$reg*.tif", TIDAL_DATA_DIR))
    if !isfile(criteria_paths[:HighTide])
        # Disk-based resample is the normal path (GDAL streams block-by-block),
        # not a fallback for OOM. The in-memory branch this replaced applied a
        # `Rasters.crop(...; to=bathy_gda2020)` after resampling `to=bathy_gda2020`,
        # which is a no-op since the result is already at that extent.
        resample_to_disk(
            Raster(hightide_fn; lazy=true), criteria_paths[:HighTide];
            to=bathy_gda2020, method=:bilinear
        )
        force_gc_cleanup()
    end

    @debug "$(now()) - Processing $(reg) - Low Tide"
    lowtide_fn = first(glob("*_lowtide_*_$reg*.tif", TIDAL_DATA_DIR))
    if !isfile(criteria_paths[:LowTide])
        # Disk-based resample is the normal path (GDAL streams block-by-block),
        # not a fallback for OOM. The in-memory branch this replaced applied a
        # `Rasters.crop(...; to=bathy_gda2020)` after resampling `to=bathy_gda2020`,
        # which is a no-op since the result is already at that extent.
        resample_to_disk(
            Raster(lowtide_fn; lazy=true), criteria_paths[:LowTide];
            to=bathy_gda2020, method=:bilinear
        )
        force_gc_cleanup()
    end

    # Calculate distance to nearest port
    @debug "$(now()) - Processing $(reg) - Ports"
    within_port_range(
        criteria_paths[:Depth],
        port_buffer,
        criteria_paths[:PortDistSlopes],
    )

    write_region_outputs(reg, criteria_paths)

    return nothing
end

"""
    prep_MPA(; parallel::Bool=false)::Nothing

Prepare data for analysis by processing UQ-GBRMPA files for each GBRMPA
management region:
1. Load the region boundary geometries.
2. Prepare inputs shared across all regions (zoning exclusion, ports, port
   buffer).
3. Process each region (see `process_region`), writing that region's
   outputs (see `write_region_outputs`).

# Arguments
- `parallel`: if `true` and additional `Distributed.jl` worker processes are
  already available (added by the caller via `addprocs`, e.g.
  `addprocs(4); prep_MPA(; parallel=true)`), regions are processed
  concurrently with `pmap`, with each worker rebuilding its own GEOS/GDAL
  state from disk (see `process_region`). Worker count is intentionally not
  chosen here — callers control it via `addprocs`. If no extra workers are
  available, falls back to the sequential loop regardless of this flag.
"""
function prep_MPA(; parallel::Bool=false)::Nothing
    # Fail fast if the region boundaries file is missing, before any of the
    # (expensive) per-region work below begins. `process_region` reloads
    # this itself, since it must also be safe to call from a remote worker.
    load_region_geometries()

    port_buffer_fn = prep_ports_and_zoning()

    # 2. Process MPA files to represent GBRMPA regions in GDA2020 projection
    # If a file already exists it is skipped (see `process_region` /
    # `write_region_outputs`).
    if parallel && nprocs() > 1
        # Make sure every worker has the function/constant definitions this
        # file provides (each worker rebuilds its own GEOS/GDAL state from
        # disk inside `process_region`, so no live handles are sent over the
        # wire - only `reg::String` and `port_buffer_fn::String`).
        @everywhere include($(@__FILE__))
        pmap(reg -> process_region(reg, port_buffer_fn), REGIONS)
    else
        if parallel
            @warn "prep_MPA(; parallel=true) requested but no extra Distributed.jl workers are available (nprocs() == 1). Add workers first, e.g. `using Distributed; addprocs(4)`, before calling `prep_MPA(; parallel=true)`. Falling back to sequential processing."
        end

        @showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
            process_region(reg, port_buffer_fn)
        end
    end

    return nothing
end
