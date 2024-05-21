"""Identify suitable locations for each region."""

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

    function assess_region(reg)
        # Load prepared bathymetry raster
        src_bathy = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"); crs=EPSG(7844), lazy=true)

        # Load prepared slope raster
        src_slope = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"); crs=EPSG(7844), lazy=true)

        # Load prepared benthic data
        src_benthic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"); crs=EPSG(7844), lazy=true)

        # Load prepared geomorphic data
        src_geomorphic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"); crs=EPSG(7844), lazy=true)

        # Load prepared wave height data
        src_waves_Hs = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"); crs=EPSG(7844), lazy=true)

        # Load prepared wave period data
        src_waves_Tp = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"); crs=EPSG(7844), lazy=true)

        # Apply filtering criteria to raster grid
        suitable_areas = read(
            (src_benthic .∈ [MPA_BENTHIC_IDS]) .&
            (-9.0 .<= src_bathy .<= -2.0) .&
            (0.0 .<= src_slope .<= 40.0) .&
            (0.0 .<= src_waves_Hs .<= 1.0) .&
            (0.0 .<= src_waves_Tp .<= 6.0)
        )

        src_bathy = nothing
        src_slope = nothing
        src_benthic = nothing
        src_waves_Hs = nothing
        src_waves_Tp = nothing

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Int16, copy(suitable_areas))
        rebuild(result_raster; missingval=0)

        # Assess flats
        suitable_flats = read(suitable_areas .& (src_geomorphic .∈ [MPA_FLAT_IDS]))

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), border=Fill(0))
        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_flats.tif")
        _write_data(fpath, res, result_raster)

        suitable_flats = nothing
        res = nothing
        GC.gc()

        # Assess slopes
        suitable_slopes = read(suitable_areas .& (src_geomorphic .∈ [MPA_SLOPE_IDS]))

        # Calculate suitability of 10x10m surroundings of each cell
        res = mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), border=Fill(0))
        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
        _write_data(fpath, res, result_raster)

        suitable_areas = nothing
        suitable_slopes = nothing
        res = nothing
        GC.gc()
    end
end

@showprogress dt = 10 desc = "Analyzing..." pmap(assess_region, REGIONS)
