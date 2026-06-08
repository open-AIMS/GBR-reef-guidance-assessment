"""
Creating a hybrid data product combining UQ-GBRMPA dataset with ACA to gap fill missing
(in-shore) areas of the GBR.

This script creates the hybridization of the two data products for benthic data.
"""

include("common.jl")

function prep_hybrid_benthic()
    management_zones = GDF.read(
        joinpath(
            CONFIG["gda2020_data"]["GDA2020_DATA_DIR"],
            "Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg"
        )
    )

    mpa_benthic_data = Raster(
        joinpath(CONFIG["mpa_data"]["MPA_DATA_DIR"], "benthic", "GBR10 GBRMP Benthic.tif");
        lazy=true,
        missingval=0
    )

    # Cache reprojected+filtered ACA polygons as GeoPackage to skip slow GeoJSON parsing
    # and reprojection on subsequent runs
    aca_cache_fn = joinpath(CONFIG["aca_data"]["ACA_DATA_DIR"], "Benthic-Map", "cache_benthic_GDA2020.gpkg")
    if isfile(aca_cache_fn)
        @debug "Loading cached ACA benthic data"
        target_polys = GDF.read(aca_cache_fn)
    else
        target_polys = GDF.read(
            joinpath(CONFIG["aca_data"]["ACA_DATA_DIR"], "Benthic-Map", "benthic.geojson")
        )

        @debug "Reprojecting ACA to target CRS"
        @time target_polys = GDF.reproject(target_polys, GI.crs(target_polys), EPSG_7844)

        @debug "Building STR tree for ACA"
        tree = STRT.STRtree(target_polys.geometry)
        reg_poly_idx = vcat(STRT.query.(Ref(tree), management_zones.SHAPE)...)
        target_polys = target_polys[reg_poly_idx, :]

        # Standardize text
        target_polys.class .= lowercase.(
            replace.(target_polys.class, " " => "_", "/" => "_")
        )

        target_polys.class_id = map(
            x -> Symbol(x) in keys(MPA_FULL_BENTHIC_IDS) ? getindex(MPA_FULL_BENTHIC_IDS, Symbol(x)) : 0,
            target_polys.class
        )

        @debug "Caching reprojected ACA benthic data"
        GDF.write(aca_cache_fn, target_polys; crs=EPSG_7844)
    end

    # Rebuild query tree
    tree = STRT.STRtree(target_polys.geometry)

    @info "Prepping hybrid benthic data"
    for reg in REGIONS
        fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_hybrid_benthic.tif")
        if isfile(fn)
            @info "Skipping $reg as file already exists..."
            continue
        end

        force_gc_cleanup()

        reg_idx = occursin.(reg[1:3], management_zones.AREA_DESCR)
        r = management_zones[reg_idx, :]

        @info "Cropping raster to management zone $reg"
        @time cropped_gbr10 = Rasters.trim(
            Rasters.mask(
                Rasters.crop(mpa_benthic_data; to=r.SHAPE);
                with=r.SHAPE
            )
        )

        # @info "Marking valid areas"
        # cropped_mz = read(cropped_mz .* Bool.(cropped_mz .∈ Ref(values(MPA_FULL_BENTHIC_IDS))))

        # Select benthic classes of interest and reproject to target CRS to ensure alignment
        @info "Ensuring reprojection is EPSG:7844"
        tmp_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic_tmp.tif")
        @time Rasters.resample(
            cropped_gbr10;
            crs=EPSG_7844,
            filename=tmp_fn
        )
        cropped_gbr10 = Raster(tmp_fn; lazy=true, missingval=0)

        reg_poly_idx = vcat(STRT.query.(Ref(tree), r.SHAPE)...)
        reg_polys = target_polys[reg_poly_idx, :]

        @info "Rasterizing ACA polygons"
        @time cropped_aca = Rasters.rasterize(
            maximum,
            reg_polys;
            to=cropped_gbr10,
            fill=:class_id,
            missingval=0
        )

        # Where GBR10 has data use it, otherwise fall back to ACA (fused, single-pass)
        # Mask ACA to management zone boundary to avoid filling bounding-box border areas
        cropped_aca = mask(cropped_aca; with=r.SHAPE)
        @info "Writing hybrid data for $reg"
        @time write_cog(fn, Int8.(ifelse.(!=(0).(cropped_gbr10), cropped_gbr10, cropped_aca)))
        rm(tmp_fn)
    end

    return nothing
end
