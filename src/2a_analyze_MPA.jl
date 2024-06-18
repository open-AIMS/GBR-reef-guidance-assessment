"""Identify suitable locations for each region."""

include("common.jl")

@everywhere begin
    using Base.Threads
    using SharedArrays

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
        split_mapwindow(f::Function, img::Raster, window::Tuple{UnitRange{Int64}, UnitRange{Int64}}, no_splits::Int = 2; border::Union{Fill, Nothing}=Fill(0))::Array

    Apply a function to a window of a raster image in parallel.
    """
    function split_mapwindow(f::Function, img::Raster, window::Tuple{UnitRange{Int64}, UnitRange{Int64}}, no_splits::Int = 2; border::Union{Fill, Nothing}=Fill(0))::Array
        shared_output_array = SharedArray{Int16}(size(img))
        split_width = Int(ceil(size(img, 2)/ no_splits))

        Threads.@threads for j in 1:no_splits  # for each col
            j_start = split_width * (j - 1) + 1

            if j == no_splits
                j_end = size(img, 2)
            else
                j_end = split_width*j
            end

            img_col = img[:, j_start:j_end]     # get column buf of img
            shared_output_array[:, j_start:j_end] = mapwindow(f, img_col, window, border) # apply function to window
        end
        output_array = Array(shared_output_array)

        shared_output_array = nothing  # Not sure if SharedArray memory cleared
        GC.gc()
        return output_array
    end

    function assess_region(reg)
        if (isfile(joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_areas.tif")))
            suitable_areas = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_areas.tif"); crs=EPSG(7844), lazy=true)

            # Need a copy of raster data type to support writing to `tif`
            result_raster = convert.(Int16, copy(suitable_areas))
            rebuild(result_raster; missingval=0)

        else
            # Load prepared rasters
            src_bathy = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_bathy.tif"); crs=EPSG(7844), lazy=true)
            src_slope = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_slope.tif"); crs=EPSG(7844), lazy=true)
            src_benthic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic.tif"); crs=EPSG(7844), lazy=true)
            src_waves_Hs = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_waves_Hs.tif"); crs=EPSG(7844), lazy=true)
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

            # write suitable_areas to .tif file
            write(joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_areas.tif"), result_raster; force=true)
        end

        src_geomorphic = Raster(joinpath(MPA_OUTPUT_DIR, "$(reg)_geomorphic.tif"); crs=EPSG(7844), lazy=true)

        # Assess flats
        suitable_flats = read(suitable_areas .& (src_geomorphic .∈ [MPA_FLAT_IDS]))
        no_splits = Int(floor(size(suitable_flats,2) / 9))  #1000 # window size? ... 10x10m

        # Calculate suitability of 10x10m surroundings of each cell
        println("Mapping $(reg[1:6]) flats... split_mapwindow($(no_splits))")
        res = split_mapwindow(prop_suitable, suitable_flats, (-4:5, -4:5), no_splits; border=Fill(0))
        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_flats.tif")
        _write_data(fpath, res, result_raster)

        suitable_flats = nothing
        res = nothing
        GC.gc()

        # Assess slopes
        suitable_slopes = read(suitable_areas .& (src_geomorphic .∈ [MPA_SLOPE_IDS]))
        no_splits = Int(floor(size(suitable_slopes,2) / 9))  #1000 # window size? ... 10x10m

        suitable_areas = nothing

        # Calculate suitability of 10x10m surroundings of each cell
        println("Mapping $(reg[1:6]) slopes... split_mapwindow($(no_splits))")
        res = split_mapwindow(prop_suitable, suitable_slopes, (-4:5, -4:5), no_splits; border=Fill(0))
        fpath = joinpath(MPA_OUTPUT_DIR, "$(reg)_suitable_slopes.tif")
        _write_data(fpath, res, result_raster)

        suitable_slopes = nothing
        res = nothing
        GC.gc()
    end
end

@time @showprogress dt = 10 desc = "Analyzing..." pmap(assess_region, REGIONS)
