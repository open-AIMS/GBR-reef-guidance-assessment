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

    function analyze_allen(reg)
        # Load required prepared raster files for analysis
        bathy_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"); crs=EPSG(7844), lazy=true)

        turbid_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"); crs=EPSG(7844), lazy=true)

        waves_Hs_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"); crs=EPSG(7844), lazy=true)

        waves_Tp_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_waves_Tp.tif"); crs=EPSG(7844), lazy=true)

        # Apply criteria to raster grid
        suitable_raster = read(
            (200 .<= bathy_rst .<= 900) .&
            (0 .<= turbid_rst .<= 52) .&
            (0.0 .<= waves_Hs_rst .<= 1.0) .&
            (0.0 .<= waves_Tp_rst .<= 6.0)
        )

        bathy_rst = nothing
        turbid_rst = nothing
        waves_Hs_rst = nothing
        waves_Tp_rst = nothing

        # Create reef-scale raster with suitable bathy, turbidity and benthic criteria
        benthic_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg).gpkg"))
        suitable_benthic = Rasters.trim(mask(suitable_raster; with=benthic_poly, boundary=:touches))

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
