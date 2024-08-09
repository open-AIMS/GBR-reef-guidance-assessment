"""
Identify suitable locations based on bathymetric, turbidity, geomorphic and benthic criteria.
Use ACA raster and polygon data to analyze suitabile areas.
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

        return Int16(round((total / length(subsection)) * 100))
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
        raster_lon = Vector{Float64}(tmp_areas.dims[1].val)
        raster_lat = Vector{Float64}(tmp_areas.dims[2].val)

        @floop for (lon_ind, lon) in enumerate(raster_lon)
            for (lat_ind, lat) in enumerate(raster_lat)
                if tmp_areas.data[lon_ind, lat_ind]

                    point = AG.createpoint()
                    AG.addpoint!(point, lon, lat)

                    pixel_dists = AG.distance.([point], Ports.geometry)
                    geom_point = gdf[argmin(pixel_dists), :geometry]
                    geom_point = (AG.getx(geom_point, 0), AG.gety(geom_point, 0))

                    dist_nearest = Distances.haversine(geom_point, (lon, lat))

                    # convert from meters to nautical miles
                    if units == "NM"
                        dist_nearest = dist_nearest / 1852
                    end

                    # convert from meters to kilometers
                    if units == "km"
                        dist_nearest = dist_nearest / 1000
                    end

                    if dist_nearest < dist
                        tmp_areas.data[lon_ind, lat_ind] = 1
                    else
                        tmp_areas.data[lon_ind, lat_ind] = 0
                    end
                end
            end
        end

        return tmp_areas
    end

    # Load QLD_ports data
    ports = GDF.read("$(PORT_DATA_DIR)/ports_QLD_merc.shp")
    ports.geometry = AG.reproject(Ports.geometry, crs(Ports[1, :geometry]), GDA2020_crs; order=:trad)

    function analyze_allen(reg)
        # Load required prepared raster files for analysis
        bathy_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"); crs=EPSG(7844), lazy=true)

        turbid_rst = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_turbid.tif"); crs=EPSG(7844), lazy=true)

        waves_Hs_rst = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"); crs=EPSG(7844), lazy=true)

        waves_Tp_rst = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"); crs=EPSG(7844), lazy=true)

        if reg == "Townsville-Whitsunday"
            rugosity_rst = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_rugosity.tif"); crs=EPSG(7844), lazy=true)

            # Apply criteria to raster grid
            suitable_raster = read(
                (200 .<= bathy_rst .<= 900) .&
                (0 .<= turbid_rst .<= 58) .&
                (0.0 .<= waves_Hs_rst .<= 1.0) .&
                (0.0 .<= waves_Tp_rst .<= 6.0) .&
                (rugosity_rst .< 6)
            )
        else
            # Apply criteria to raster grid
            suitable_raster = read(
                (200 .<= bathy_rst .<= 900) .&
                (0 .<= turbid_rst .<= 58) .&
                (0.0 .<= waves_Hs_rst .<= 1.0) .&
                (0.0 .<= waves_Tp_rst .<= 6.0)
            )
        end

        bathy_rst = nothing
        turbid_rst = nothing
        waves_Hs_rst = nothing
        waves_Tp_rst = nothing
        rugosity_rst = nothing

        # Create reef-scale raster with suitable bathy, turbidity and benthic criteria
        benthic_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg).gpkg"))
        suitable_benthic = Rasters.trim(mask(suitable_raster; with=benthic_poly, boundary=:touches))

        suitable_benthic = filter_distances(suitable_benthic, ports, 200; units="NM")

        # Filter out cells occurring
        GBRMPA_zone_exclusion = GDF.read(joinpath(MPA_OUTPUT_DIR, "GBRMPA_zone_exclusion.gpkg"))
        suitable_benthic = Rasters.mask(suitable_benthic; with=GBRMPA_zone_exclusion, invert=true, boundary=:touches)

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_benthic))
        rebuild(result_raster; missingval=0)

        # Asses flats
        target_flats_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg).gpkg"))
        suitable_flats = Rasters.mask(suitable_benthic; with=target_flats_poly, boundary=:touches)

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), border=Fill(0))
        fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_flats.tif")
        _write_data(fpath, res, result_raster)

        suitable_flats = nothing
        res = nothing
        GC.gc()

        # Assess slopes
        target_slopes_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg).gpkg"))
        suitable_slopes = Rasters.mask(suitable_benthic; with=target_slopes_poly, boundary=:touches)

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), border=Fill(0))
        fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
        _write_data(fpath, res, result_raster)

        suitable_slopes = nothing
        res = nothing
        result_raster = nothing
        suitable_benthic = nothing
        GC.gc()
    end
end

@showprogress dt = 10 desc = "Analyzing..." pmap(analyze_allen, REGIONS)
