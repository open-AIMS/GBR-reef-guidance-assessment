"""Identify suitable locations for each region."""

include("common.jl")


@everywhere begin
    """
        prop_suitable(subsection::AbstractMatrix)::Float32

    Calculate the the proportion of the subsection that is suitable for deployments.
    Subsection is the surrounding hectare centred on each cell of a raster.
    """
    function prop_suitable(subsection::AbstractMatrix)::Float32
        total = sum(subsection)
        if total == 0.0
            return 0.0
        end

        return Float32((total / length(subsection)))
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

    function assess_region(reg)
        # Load bathymetry raster
        src_bathy_path = joinpath(OUTPUT_DIR, "$(reg)_bathy.tif")
        src_bathy = Raster(src_bathy_path; crs=EPSG(7844), lazy=true)

        # Load slope raster
        src_slope_path = joinpath(OUTPUT_DIR, "$(reg)_slope.tif")
        src_slope = Raster(src_slope_path; crs=EPSG(7844), lazy=true)

        # Load pre-prepared benthic data
        src_benthic_path = joinpath(OUTPUT_DIR, "$(reg)_benthic.tif")
        src_benthic = Raster(src_benthic_path; crs=EPSG(7844), lazy=true)

        # Load pre-prepared geomorphic data
        src_geomorphic_path = joinpath(OUTPUT_DIR, "$(reg)_geomorphic.tif")
        src_geomorphic = Raster(src_geomorphic_path; crs=EPSG(7844), lazy=true)

        # Load pre-prepared wave height data
        src_waves_Hs_path = joinpath(OUTPUT_DIR, "$(reg)_waves_Hs.tif")
        src_waves_Hs = Raster(src_waves_Hs_path; crs=EPSG(7844), lazy=true)

        # Load pre-prepared wave period data
        src_waves_Tp_path = joinpath(OUTPUT_DIR, "$(reg)_waves_Tp.tif")
        src_waves_Tp = Raster(src_waves_Tp_path; crs=EPSG(7844), lazy=true)

        # Source image is of 10m^2 pixels
        # A hectare is 100x100 meters, so we calculate the proportional area of each hectare
        # that meet criteria of (-9 <= depth <= -3m, slope < 40 deg, habitat is Rock or Coral/Algae,
        # 90th percentile of standing wave height is below 1m, and wave period is less than 6 sec).

        suitable_areas = read(
            (src_benthic .∈ [BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0) .&
            (0.0 .<= src_waves_Hs .<= 1.0) .&
            (0.0 .<= src_waves_Tp .<= 6.0)
        )

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Float32, copy(suitable_areas))
        rebuild(result_raster; missingval=0)

        # Assess flats
        suitable_flats = read((suitable_areas) .& (src_geomorphic .∈ [FLAT_IDS]))

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_flats.tif")
        _write_data(fpath, res, result_raster)

        suitable_flats = nothing
        res = nothing
        GC.gc()

        # Assess slopes
        suitable_slopes = read((suitable_areas) .& (src_geomorphic .∈ [SLOPE_IDS]))

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
        _write_data(fpath, res, result_raster)

        suitable_slopes = nothing
        res = nothing
        GC.gc()
    end
end

@showprogress dt = 10 desc = "Analyzing..." map(assess_region, REGIONS)
