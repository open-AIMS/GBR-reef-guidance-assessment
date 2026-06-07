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

# This processing step requires all the memory available so we actually remove the workers
# that were just set up.
# This can be commented out if working on a system with large enough amounts of memory
rmprocs(workers()[2:end]...)
GC.gc()

using SparseArrays, NamedTupleTools
using ExtendableSparse

# Loading regions_4326 for cropping of vector and raster data.
regions_4326 = GDF.read(REGION_PATH_4326)

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

# 1c. Create buffer around ports (200 NM said to be max range considered)
if !isfile(joinpath(MPA_OUTPUT_DIR, "ports_buffer.gpkg"))
    port_locs = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_GDA2020.gpkg"))

    port_buffer = port_buffer_mask(port_locs, 200.0, unit="NM")
    port_buffer = DataFrame(; Name="ports_buffer", geometry=port_buffer)
    GDF.write(
        joinpath(MPA_OUTPUT_DIR, "port_buffer.gpkg"),
        port_buffer;
        crs=EPSG_7844
    )
end


# 2. Process MPA files to represent GBRMPA regions in GDA2020 projection

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

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
            -9999.0,
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
            -9999.0,
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
            Raster(raw_rugosity_fn; crs=REGION_CRS_UTM[reg], mappedcrs=EPSG_4326),
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
        -9999.0;
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
        -9999.0;
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
        -9999.0;
        method=:bilinear
    )

    # Tidal data is already separated into management regions
    # so we only need to reproject into consistent datum
    @debug "$(now()) - Processing $(reg) - High Tide"
    hightide_fn = first(glob("*_hightide_*_$reg*.tif", TIDAL_DATA_DIR))
    if !isfile(criteria_paths[:HighTide])
        high_tide = resample(Raster(hightide_fn); to=bathy_gda2020, method=:bilinear)
        Rasters.write(
            criteria_paths[:HighTide],
            Rasters.crop(high_tide; to=bathy_gda2020)
        )

        high_tide = nothing
        force_gc_cleanup()
    end

    @debug "$(now()) - Processing $(reg) - Low Tide"
    lowtide_fn = first(glob("*_lowtide_*_$reg*.tif", TIDAL_DATA_DIR))
    if !isfile(criteria_paths[:LowTide])
        low_tide = resample(Raster(lowtide_fn); to=bathy_gda2020, method=:bilinear)
        Rasters.write(
            criteria_paths[:LowTide],
            Rasters.crop(low_tide; to=bathy_gda2020)
        )

        low_tide = nothing
        force_gc_cleanup()
    end

    # Calculate distance to nearest port
    @debug "$(now()) - Processing $(reg) - Ports"
    within_port_range(
        criteria_paths[:Depth],
        port_buffer,
        criteria_paths[:PortDistSlopes],
    )

    @debug "$(now()) - Processing $(reg) - Valid area"
    # Find locations containing valid data
    valid_slopes_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif")
    if !isfile(valid_slopes_fn)
        write_valid_locs(
            criteria_paths,
            [values(MPA_BENTHIC_IDS)...],
            [values(MPA_SLOPE_IDS)...],
            7, (3, 3), 70, (9, 9),
            valid_slopes_fn
        )

        resize_to_valid_area(
            criteria_paths,
            valid_slopes_fn
        )
    end

    # Create lookup tables to support fast querying
    @debug "$(now()) - Processing $(reg) - Lookup table"
    slopes_lookup_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes_lookup.parq")
    valid_lookup(
        criteria_paths,
        valid_slopes_fn,
        slopes_lookup_fn
    )

    # valid_flats_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats.tif")
    # write_valid_locs(
    #     criteria_paths,
    #     MPA_BENTHIC_IDS,
    #     MPA_FLAT_IDS,
    #     7, (3, 3), 70, (9, 9),
    #     valid_flats_fn,
    #     reg
    # )

    # flats_lookup_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_flats_lookup.parq")
    # valid_lookup(
    #     criteria_paths,
    #     valid_flats_fn,
    #     flats_lookup_fn
    # )
end
