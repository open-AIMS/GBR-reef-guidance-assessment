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
        crs=EPSG(7844)
    )
end

# 2. Process MPA files to represent GBRMPA regions in GDA2020 projection

"""
    set_consistent_missingval(raster, val)

Replace value used to indicate no data, and return a Raster type with this value set.

## Note
`raster` should cover the smallest possible extent to reduce processing time.

# Arguments
- `raster` : data to reconfigure
- `val` : Value to use to indicate no data

# Returns
Raster
"""
function set_consistent_missingval(raster, val)
    replace_missing!(raster, val)
    return Raster(raster; missingval=val)
end

# If a file already exists it is skipped
@showprogress dt = 10 "Prepping benthic/geomorphic/wave data..." for reg in REGIONS
    reg_idx_4326 = occursin.(reg[1:3], regions_4326.AREA_DESCR)

    base_bathy_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif")
    if !isfile(base_bathy_fn)
        target_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        target_bathy = Raster(target_bathy_path; mappedcrs=EPSG_4326)
        target_bathy = set_consistent_missingval(target_bathy, -9999.0)

        write(
            base_bathy_fn,
            resample(
                target_bathy;
                crs=GDA2020_crs
            )
        )

        target_bathy = nothing
        force_gc_cleanup()
    end

    # Load bathymetry data to provide corresponding spatial extent
    bathy_gda2020 = Raster(base_bathy_fn; crs=EPSG_7844, lazy=true)

    slope_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif")
    if !isfile(slope_fn)
        target_slope_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        target_slope = Raster(target_slope_path, mappedcrs=EPSG_4326)
        target_slope = set_consistent_missingval(target_slope, -9999.0)

        write(
            slope_fn,
            resample(
                target_slope;
                crs=GDA2020_crs
            )
        )

        target_slope = nothing
        force_gc_cleanup()
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

        write(
            benthic_fn,
            resample(
                target_benthic;
                to=bathy_gda2020
            )
        )

        gbr_benthic = nothing
        target_benthic = nothing
        force_gc_cleanup()
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

        write(
            geomorphic_fn,
            resample(
                target_geomorphic;
                to=bathy_gda2020
            )
        )

        target_geomorphic = nothing
        force_gc_cleanup()
    end

    target_waves_Hs_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif")
    if !isfile(target_waves_Hs_fn)
        waves_Hs_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Hs", reg)))

        # Have to load netCDF data into memory to allow missing value replacement
        waves_Hs = Raster(
            waves_Hs_path,
            key=:Hs90,
            crs=GDA2020_crs,
            mappedcrs=EPSG_4326
        )

        # Manually set -infinite missing data value to exact value
        # This is necessary as the netCDF was provided without a set `no data` value
        waves_Hs.data[waves_Hs.data .< -9999.0] .= -9999.0

        # Reset data type from Union{Missing,Float32} -> Float32
        # and also set missing value
        waves_Hs = Raster(
            waves_Hs;
            data=Matrix{Float32}(waves_Hs.data),
            missingval=-9999.0
        )

        # This NetCDF is in reverse orientation, so we have to flip it back.
        # (remember that the Y dimension is the columns, which corresponds to the longitude)
        write(
            target_waves_Hs_fn,
            reverse(
                resample(
                    waves_Hs;
                    to=bathy_gda2020
                );
                dims=Y
            )
        )

        waves_Hs = nothing
        force_gc_cleanup()
    end

    waves_Tp_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif")
    if !isfile(waves_Tp_fn)
        target_waves_Tp_path = first(glob("*.nc", joinpath(WAVE_DATA_DIR, "Tp", reg)))

        # Have to load netCDF data into memory to allow missing value replacement
        waves_Tp = Raster(
            target_waves_Tp_path,
            key=:Tp90,
            crs=GDA2020_crs,
            mappedcrs=EPSG_4326
        )
        waves_Tp.data[waves_Tp.data .< -9999.0] .= -9999.0

        # Reset data type from Union{Missing,Float32} -> Float32
        # and also set missing value
        waves_Tp = Raster(
            waves_Tp;
            data=Matrix{Float32}(waves_Tp.data),
            missingval=-9999.0
        )

        # This NetCDF is in reverse orientation, so we have to flip it back.
        # (remember that the Y dimension is the columns, which corresponds to the longitude)
        write(
            waves_Tp_fn,
            reverse(
                resample(
                    waves_Tp;
                    to=bathy_gda2020
                );
                dims=Y
            )
        )

        waves_Tp = nothing
        force_gc_cleanup()
    end

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

        write(
            turbid_fn,
            resample(
                target_turbid;
                to=bathy_gda2020
            )
        )

        target_turbid = nothing
        force_gc_cleanup()
    end

    if reg == "Townsville-Whitsunday"
        rugosity_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif")
        if !isfile(rugosity_fn)
            tsv_rugosity_path = joinpath(RUG_DATA_DIR, "std25_Rugosity_Townsville-Whitsunday.tif")
            tsv_rugosity = Raster(tsv_rugosity_path, mappedcrs=EPSG_4326)
            tsv_rugosity = set_consistent_missingval(tsv_rugosity, -9999.0)

            write(
                rugosity_fn,
                resample(
                    tsv_rugosity;
                    to=bathy_gda2020
                )
            )

            tsv_rugosity = nothing
            force_gc_cleanup()
        end
    end
end
