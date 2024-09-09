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
include("geom_handlers/raster_processing.jl")
include("geom_handlers/lookup_processing.jl")

# This processing step requires all the memory available so we actually remove the workers
# that were just set up.
# This can be commented out if working on a system with large enough amounts of memory
rmprocs(workers()[2:end]...)
GC.gc()

using SparseArrays, NamedTupleTools

# Loading regions_4326 for cropping of vector and raster data.
regions_4326 = GDF.read(REGION_PATH_4326)

# 1. Process GBRMPA zoning geopackage to select only zones for site exclusion
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

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    # Create NamedTuple to hold all output file paths.
    criteria_paths = (
        Depth = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"),
        Benthic = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"),
        Geomorphic = joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"),
        Slope = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"),
        Turbidity = joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"),
        WavesHs = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"),
        WavesTp = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"),
        Rugosity = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"),
        PortDistSlopes = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_slopes.tif"),
        PortDistFlats = joinpath(MPA_OUTPUT_DIR, "$(reg)_port_distance_flats.tif")
    )
    if reg != "Townsville-Whitsunday"
        criteria_paths = NamedTupleTools.delete(criteria_paths, :Rugosity)
    end

    # Process bathymetry, slope and rugosity UTM raster files
    raw_bathy_fn = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    process_UTM_raster(raw_bathy_fn, criteria_paths[:Depth], EPSG_7844, -9999.0, reg)

    raw_slope_fn = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
    process_UTM_raster(raw_slope_fn, criteria_paths[:Slope], EPSG_7844, -9999.0, reg)

    if reg == "Townsville-Whitsunday"
        raw_rugosity_fn = joinpath(RUG_DATA_DIR, "std25_Rugosity_Townsville-Whitsunday.tif")
        process_UTM_raster(raw_rugosity_fn, criteria_paths[:Rugosity], EPSG_7844, -9999.0, reg)
    end

    # Process GBR-wide raster data
    # Load bathymetry data to provide corresponding spatial extent
    bathy_gda2020 = Raster(criteria_paths[:Depth]; crs=EPSG_7844, lazy=true)

    raw_benthic_fn = "$(MPA_DATA_DIR)/benthic/GBR10 GBRMP Benthic.tif"
    target_benthic = trim_extent_region(
        raw_benthic_fn,
        EPSG_4326,
        regions_4326[reg_idx_4326, :geometry],
        criteria_paths[:Benthic]
    )
    resample_and_write(target_benthic, bathy_gda2020, criteria_paths[:Benthic])
    target_benthic = nothing
    force_gc_cleanup()

    raw_geomorphic_fn = "$(MPA_DATA_DIR)/geomorphic/GBR10 GBRMP Geomorphic.tif"
    target_geomorphic = trim_extent_region(
        raw_geomorphic_fn,
        EPSG_4326,
        regions_4326[reg_idx_4326, :geometry],
        criteria_paths[:Geomorphic]
    )
    resample_and_write(target_geomorphic, bathy_gda2020, criteria_paths[:Geomorphic])
    target_geomorphic = nothing
    force_gc_cleanup()

    raw_turbid_fn = "$(ACA_DATA_DIR)/Turbidity-Q3-2023/turbidity-quarterly_0.tif"
    target_turbid = trim_extent_region(
        raw_turbid_fn,
        EPSG_4326,
        regions_4326[reg_idx_4326, :geometry],
        criteria_paths[:Turbidity]
    )
    resample_and_write(target_turbid, bathy_gda2020, criteria_paths[:Turbidity])
    target_turbid = nothing
    force_gc_cleanup()

    # Process wave raster data
    # Use bathy dataset as a template for writing netCDF data to geotiff
    src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
    rst_template = Raster(src_bathy_path, crs=REGION_CRS_UTM[reg], mappedcrs=EPSG(4326), lazy=true)

    waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))
    process_wave_data(waves_Hs_path, criteria_paths[:WavesHs], :Hs90, rst_template, bathy_gda2020, -9999.0)

    waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Tp", reg)))
    process_wave_data(waves_Tp_path, criteria_paths[:WavesTp], :Tp90, rst_template, bathy_gda2020, -9999.0)

    # Find locations containing valid data
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

    distance_raster(
        valid_slopes_fn,
        port_buffer,
        port_points,
        -9999.0,
        criteria_paths[:PortDistSlopes],
        "NM"
    )

    distance_raster(
        valid_flats_fn,
        port_buffer,
        port_points,
        -9999.0,
        criteria_paths[:PortDistFlats],
        "NM"
    )

    # Create lookup tables to support fast querying
    slopes_lookup_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes_lookup.parq")
    valid_lookup(
        NamedTupleTools.delete(criteria_paths, :PortDistFlats),
        valid_slopes_fn,
        slopes_lookup_fn
    )
    flats_lookup_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats_lookup.parq")
    valid_lookup(
        NamedTupleTools.delete(criteria_paths, :PortDistSlopes),
        valid_flats_fn,
        flats_lookup_fn
    )
end
