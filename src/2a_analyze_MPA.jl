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
        src_bathy_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "bathy", reg)))
        src_bathy = Raster(src_bathy_path, mappedcrs=EPSG(4326), lazy=true)

        # Load slope raster
        src_slope_path = first(glob("*.tif", joinpath(MPA_DATA_DIR, "slope", reg)))
        src_slope = Raster(src_slope_path, mappedcrs=EPSG(4326), lazy=true)

        # Load pre-prepared benthic data
        src_benthic_path = joinpath(OUTPUT_DIR, "$(reg)_benthic.tif")
        src_benthic = Raster(src_benthic_path, lazy=true)

        # Load pre-prepared geomorphic data
        src_geomorphic_path = joinpath(OUTPUT_DIR, "$(reg)_geomorphic.tif")
        src_geomorphic = Raster(src_geomorphic_path, lazy=true)

        # Load pre-prepared wave data
        src_waves_path = joinpath(OUTPUT_DIR, "$(reg)_waves.tif")
        src_waves = Raster(src_waves_path, lazy=true, crs=crs(src_bathy))

        # Source image is of 10m^2 pixels
        # A hectare is 100x100 meters, so we calculate the proportional area of each hectare
        # that meet criteria of (-9 <= depth <= -3, slope < 40, habitat is Rock or
        # Coral/Algae and 90th percentile of standing wave height is below 1m).

        suitable_areas = read(
            (src_benthic .∈ [BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0) .&
            (0.0 .<= src_waves .<= 1.0)
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
