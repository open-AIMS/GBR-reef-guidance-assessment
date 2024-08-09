"""
Identify suitable locations based on bathymetric, slope, geomorphic, benthic and wave criteria.
Use MPA raster data to analyze suitabile areas.
Output raster files with suitability proportions.
"""

include("common.jl")

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

    """
        filter_distances(
            target_rast::Raster,
            gdf::DataFrame,
            dist_nm
        )::Raster

    Exclude pixels in target_rast that are beyond `dist_nm` (nautical miles) from a geometry
    in `gdf`. Target_rast and gdf should be in the same CRS (EPSG:7844 / GDA2020 for GBR-reef-guidance-assessment).

    # Arguments
    - `target_rast` : Raster of suitable pixels (Bool) to filter pixels from.
    - `gdf` : DataFrame with `geometry` column that contains vector objects of interest.
    - `dist_nm` : Filtering distance from geometry object in nautical miles.

    # Returns
    - `tmp_areas` : Raster of filtered pixels containing only pixels within target distance
    from a geometry centroid.
    """
    function filter_distances(target_rast::Raster, gdf::DataFrame, dist; units::String="NM")::Raster
        tmp_areas = copy(target_rast)

        # First dimension is the rows (latitude)
        # Second dimension is the cols (longitude)
        raster_lat = Vector{Float64}(tmp_areas.dims[1].val)
        raster_lon = Vector{Float64}(tmp_areas.dims[2].val)

        @floop for row_col in findall(tmp_areas)
            (lat_ind, lon_ind) = Tuple(row_col)
            point = AG.createpoint()

            lon = raster_lon[lon_ind]
            lat = raster_lat[lat_ind]
            AG.addpoint!(point, lon, lat)

            pixel_dists = AG.distance.([point], port_locs.geometry)
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

            tmp_areas.data[lon_ind, lat_ind] = dist_nearest < dist ? 1 : 0
        end

        return tmp_areas
    end

    # Load QLD_ports data
    Ports = GDF.read("$(PORT_DATA_DIR)/ports_QLD_merc.shp")
    Ports.geometry = AG.reproject(Ports.geometry, crs(Ports[1, :geometry]), GDA2020_crs; order=:trad)

    function assess_region(reg)
        # Load required prepared raster files for analysis
        src_bathy = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"); crs=EPSG(7844), lazy=true)

        src_slope = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"); crs=EPSG(7844), lazy=true)

        src_benthic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"); crs=EPSG(7844), lazy=true)

        src_geomorphic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"); crs=EPSG(7844), lazy=true)

        src_waves_Hs = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"); crs=EPSG(7844), lazy=true)

        src_waves_Tp = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"); crs=EPSG(7844), lazy=true)

        src_turbid = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"); crs=EPSG(7844), lazy=true)

        if reg == "Townsville-Whitsunday"
            src_rugosity = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"); crs=EPSG(7844), lazy=true)

            suitable_areas = read(
                (src_benthic .∈ [MPA_BENTHIC_IDS]) .&
                (-9.0 .<= src_bathy .<= -2.0) .&
                (0.0 .<= src_slope .<= 40.0) .&
                (0.0 .<= src_waves_Hs .<= 1.0) .&
                (0.0 .<= src_waves_Tp .<= 6.0) .&
                (src_turbid .<= 58) .&
                (src_rugosity .< 6)
            )
        else
            # Apply filtering criteria to raster grid
            suitable_areas = read(
                (src_benthic .∈ [MPA_BENTHIC_IDS]) .&
                (-9.0 .<= src_bathy .<= -2.0) .&
                (0.0 .<= src_slope .<= 40.0) .&
                (0.0 .<= src_waves_Hs .<= 1.0) .&
                (0.0 .<= src_waves_Tp .<= 6.0) .&
                (src_turbid .<= 58)
            )
        end

        src_bathy = nothing
        src_slope = nothing
        src_benthic = nothing
        src_waves_Hs = nothing
        src_waves_Tp = nothing
        src_turbid = nothing
        src_rugosity = nothing

        # Filter out cells over 200NM from the nearest port
        suitable_areas = filter_distances(suitable_areas, Ports, 200; units="NM")

        # Filter out cells occurring
        GBRMPA_zone_exclusion = GDF.read(joinpath(MPA_OUTPUT_DIR, "GBRMPA_zone_exclusion.gpkg"))
        suitable_areas = Rasters.mask(suitable_areas; with=GBRMPA_zone_exclusion, invert=true, boundary=:touches)

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_areas))
        rebuild(result_raster; missingval=0)

        # Assess flats
        suitable_flats = read(suitable_areas .& (src_geomorphic .∈ [MPA_FLAT_IDS]))

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
        GC.gc()

        # Assess slopes
        suitable_slopes = read(suitable_areas .& (src_geomorphic .∈ [MPA_SLOPE_IDS]))

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
        GC.gc()
    end
end

@showprogress dt = 10 desc = "Analyzing..." pmap(assess_region, REGIONS)
