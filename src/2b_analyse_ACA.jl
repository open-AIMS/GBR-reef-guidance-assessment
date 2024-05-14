"""
Identify suitable locations based on bathymetric, turbidity, geomorphic and benthic criteria.
Use raster and polygon data to analyze suitabile areas.
Output raster files with suitability proportions.
"""

include("common.jl")

# # Create logger to hide Rasters.mask() info output - mask() progress bars still displayed
# no_info = Logging.ConsoleLogger(stderr, Logging.Warn)

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

    function analyze_allen(reg)
        bathy_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"), lazy=true)
        turbid_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"), lazy=true)

        # # Source image is of 10m^2 pixels
        # # A hectare is 100x100 meters, so we calculate the proportional area of each hectare
        # that meets criteria (200 <= depth <= 900 cm, and turbidty is no more than LOW (5.2 FNU))

        suitable_raster = read(
            (200.0 .<= bathy_rst .<= 900.0) .&
            (turbid_rst .<= 52)
        )

        # Create reef-scale raster with suitable bathy, turbidity and benthic criteria
        benthic_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_benthic_$(reg).gpkg"))
        suitable_benthic = Rasters.trim(mask(suitable_raster; with=benthic_poly, boundary=:touches))

        # Need a copy of raster data type to support writing to `tif`
        result_raster = convert.(Float32, copy(suitable_benthic))
        rebuild(result_raster; missingval=0)

        # Flats
        target_flats_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_target_flats_$(reg).gpkg"))
        suitable_flats = Rasters.mask(suitable_benthic; with=target_flats_poly, boundary=:touches)

        res = mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray #needs to be changed to calculate suitability for hectare not 37m^2
        fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_flats.tif")
        _write_data(fpath, res, result_raster)

        res = nothing
        GC.gc()

        # Slopes
        target_slopes_poly = GDF.read(joinpath(ACA_OUTPUT_DIR, "aca_target_slopes_$(reg).gpkg"))
        suitable_slopes = Rasters.mask(suitable_benthic; with=target_slopes_poly, boundary=:touches)

        res = mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
        fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
        _write_data(fpath, res, result_raster)

        suitable_slopes = nothing
        res = nothing
        result_raster = nothing
        suitable_benthic = nothing
        GC.gc()
    end
end

@showprogress dt = 10 desc = "Analysing..." pmap(analyze_allen, REGIONS)
