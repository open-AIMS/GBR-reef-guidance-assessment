"""
Identify suitable locations based on bathymetric, slope, geomorphic, benthic and wave criteria.
Use MPA raster data to analyze suitabile areas.
Output raster files with suitability proportions.
"""

include("common.jl")
include("geom_handlers/raster_processing.jl")
include("geom_handlers/geom_ops.jl")

@everywhere begin
    """
        prop_suitable(subsection::AbstractMatrix)::Float32

    Calculate the the proportion of the subsection that is suitable for deployments.
    Subsection is the surrounding hectare centred on each cell of a raster.
    """
    function prop_suitable(subsection::AbstractMatrix)::Int16
        total = sum(subsection)
        if total == 0.0
            return 0.0
        end

        return Int16(round((total / length(subsection))*100))
    end

    function _write_data(fpath::String, data, cache)::Nothing
        if !isfile(fpath)
            if occursin("grouped_", fpath)
                cache .= label_components(data)
            else
                cache .= data
            end

            write(fpath, cache; force=true)
        end

        return nothing
    end

    # Load QLD_ports buffer data
    port_buffer = GDF.read(joinpath(MPA_OUTPUT_DIR, "ports_buffer.gpkg"))[:,:geometry]

    function assess_region(reg, port_buffer)
        # Load required prepared raster files for analysis
        src_bathy = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"); crs=EPSG_7844)
        bathy_crit = (-9.0 .<= src_bathy .<= -2.0)
        src_bathy = nothing

        src_slope = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"); crs=EPSG_7844)
        slope_crit = (0.0 .<= src_slope .<= 40.0)
        src_slope = nothing

        src_benthic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"); crs=EPSG_7844)
        benthic_crit = (src_benthic .∈ [MPA_BENTHIC_IDS])
        src_benthic = nothing

        src_geomorphic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"); crs=EPSG_7844)
        geomorphic_flat_crit = src_geomorphic .∈ [MPA_FLAT_IDS]
        geomorphic_slope_crit = src_geomorphic .∈ [MPA_SLOPE_IDS]
        src_geomorphic = nothing

        src_waves_Hs = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"); crs=EPSG_7844)
        Hs_waves_crit = (0.0 .<= src_waves_Hs .<= 1.0)
        src_waves_Hs = nothing

        src_waves_Tp = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"); crs=EPSG_7844)
        Tp_waves_crit = (0.0 .<= src_waves_Tp .<= 6.0)
        src_waves_Tp = nothing

        src_turbid = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"); crs=EPSG_7844)
        turbid_crit = (src_turbid .<= 58)
        src_turbid = nothing

        # Apply filtering criteria to raster grid
        suitable_areas = (
            benthic_crit .&
            bathy_crit .&
            slope_crit .&
            Hs_waves_crit .&
            Tp_waves_crit .&
            turbid_crit
        )

        if reg == "Townsville-Whitsunday"
            src_rugosity = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"); crs=EPSG_7844)
            suitable_areas .= suitable_areas .& (src_rugosity .< 6)
            src_rugosity = nothing
        end

        # Filter out cells over 200NM from the nearest port
        bathy_crit = filter_distances(bathy_crit, port_buffer)

        # Filter out cells occurring in preservation zones
        GBRMPA_zone_exclusion = GDF.read(joinpath(MPA_OUTPUT_DIR, "GBRMPA_preservation_zone_exclusion.gpkg"))
        region_extent = GI.extent(suitable_areas)
        rst_extent = [
            GI.Polygon(
                create_poly(create_bbox(region_extent.X, region_extent.Y), EPSG(7844))
            )
        ]
        in_region = GO.within.(GBRMPA_zone_exclusion.geometry, rst_extent)
        GBRMPA_zone_exclusion = GBRMPA_zone_exclusion[in_region, :]

        for polygon in GBRMPA_zone_exclusion.geometry
            suitable_areas = mask(suitable_areas; with=polygon, invert=true, boundary=:touches)
        end

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_areas))
        rebuild(result_raster; missingval=0)

        # Assess flats
        suitable_flats = suitable_areas .& geomorphic_flat_crit

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), border=Fill(0))
        if reg == "Townsville-Whitsunday"
            fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_flats_rugosity.tif")
        else
            fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_flats.tif")
        end
        _write_data(fpath, res, result_raster)

        suitable_flats = nothing
        res = nothing
        force_gc_cleanup()

        # Assess slopes
        suitable_slopes = suitable_areas .& geomorphic_slope_crit

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), border=Fill(0))
        if reg == "Townsville-Whitsunday"
            fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes_rugosity.tif")
        else
            fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
        end
        _write_data(fpath, res, result_raster)

        suitable_areas = nothing
        suitable_slopes = nothing
        res = nothing
        force_gc_cleanup()
    end
end

@showprogress dt = 10 desc = "Analyzing..." map(x -> assess_region(x, port_buffer), REGIONS)
