"""
Identify suitable locations based on bathymetric, turbidity, geomorphic and benthic criteria.
Use raster and polygon data to analyze suitabile areas.
Output raster files with suitability proportions.
"""

using Logging

include("common.jl")

# Create logger to hide Rasters.mask() info output - mask() progress bars still displayed
no_info = Logging.ConsoleLogger(stderr, Logging.Warn)

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
end

function analyse_allen(reg)
    bathy_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_bathy.tif"), lazy=true)
    turbid_rst = Raster(joinpath(ACA_OUTPUT_DIR, "$(reg)_turbid.tif"), lazy=true)

    # # Source image is of 10m^2 pixels
    # # A hectare is 100x100 meters, so we calculate the proportional area fo each hectare
    # that meets criteria (200 <= depth <= 900 cm, and turbidty is no more than LOW (5.2 FNU))

    suitable_raster = read(
        (200.0 .<= bathy_rst .<= 900.0)
        .&
        (turbid_rst .<= 52)
    )

    # Create reef-scale raster with suitable bathy, turbidity, geomorphic and benthic criteria
    benthic_morphic_poly = GDF.read(first(glob("*.gpkg", ACA_OUTPUT_DIR)))
    Logging.with_logger(no_info) do
        suitable_reefs = Rasters.trim(mask(suitable_raster; with=benthic_morphic_poly))
    end

    benthic_morphic_poly = nothing
    bathy_rst = nothing
    turbid_rst = nothing
    GC.gc()

    # Find suitable benthic areas
    benthic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Benthic-Map", "benthic.geojson"))

    target_benthic = benthic_poly[benthic_poly.class.∈[ACA_BENTHIC_IDS], :]
    Logging.with_logger(no_info) do
        suitable_areas = Rasters.trim(mask(suitable_reefs; with=target_benthic))
    end

    # Need a copy of raster data type to support writing to `tif`
    result_raster = convert.(Float32, copy(suitable_areas))
    rebuild(result_raster; missingval=0)

    suitable_reefs = nothing
    GC.gc()

    # Apply geomorphic classifications
    geomorphic_poly = GDF.read(joinpath(ACA_DATA_DIR, "Geomorphic-Map", "geomorphic.geojson"))

    # flats
    target_flats = geomorphic_poly[geomorphic_poly.class.∈[ACA_FLAT_IDS], :]
    Logging.with_logger(no_info) do
        suitable_flats = Rasters.mask(suitable_areas; with=target_flats)
    end

    res = mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), border=Fill(0)) .|> Gray #needs to be changed to calculate suitability for hectare not 37m^2
    fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_flats.tif")
    _write_data(fpath, res, result_raster)

    suitable_flats = nothing
    res = nothing
    GC.gc()

    # slopes
    target_slopes = geomorphic_poly[geomorphic_poly.class.∈[ACA_SLOPE_IDS], :]
    Logging.with_logger(no_info) do
        suitable_slopes = Rasters.mask(suitable_areas; with=target_slopes)
    end

    res = mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), border=Fill(0)) .|> Gray
    fpath = joinpath(ACA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
    _write_data(fpath, res, result_raster)

    suitable_slopes = nothing
    res = nothing
    result_raster = nothing
    suitable_areas = nothing
    GC.gc()
end

@showprogress dt = 10 desc = "Analysing..." map(analyse_allen, REGIONS)
