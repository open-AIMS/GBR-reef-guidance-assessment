"""
    prep_hybrid_layer(
        layer_name::String,
        mpa_source_path::String,
        aca_subdir::String,
        id_lookup::NamedTuple
    )::Nothing

Create a hybrid raster product per region, combining the UQ-GBRMPA dataset with ACA to
gap-fill missing (in-shore) areas of the GBR. Shared implementation for both the benthic
and geomorphic hybrid layers, which differ only in source raster path, ACA geojson
subdirectory, and class ID lookup table.

# Arguments
- `layer_name` : "benthic" or "geomorphic" - used in output/cache file names and log messages.
- `mpa_source_path` : Path to the GBR-wide UQ-GBRMPA raster, relative to `MPA_DATA_DIR`.
- `aca_subdir` : ACA data subdirectory holding the geojson/cache for this layer, relative to `ACA_DATA_DIR`.
- `id_lookup` : NamedTuple mapping class name (Symbol) to raster ID, used to build `class_id`.
"""
function prep_hybrid_layer(
    layer_name::String,
    mpa_source_path::String,
    aca_subdir::String,
    id_lookup::NamedTuple
)::Nothing
    management_zones = GDF.read(
        joinpath(
            CONFIG["gda2020_data"]["GDA2020_DATA_DIR"],
            "Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg"
        )
    )

    mpa_data = Raster(
        joinpath(CONFIG["mpa_data"]["MPA_DATA_DIR"], mpa_source_path);
        lazy=true,
        missingval=0
    )

    # Cache reprojected+filtered ACA polygons as GeoPackage to skip slow GeoJSON parsing
    # and reprojection on subsequent runs
    aca_cache_fn = joinpath(
        CONFIG["aca_data"]["ACA_DATA_DIR"], aca_subdir, "cache_$(layer_name)_GDA2020.gpkg"
    )
    if isfile(aca_cache_fn)
        @debug "Loading cached ACA $(layer_name) data"
        target_polys = GDF.read(aca_cache_fn)
    else
        target_polys = GDF.read(
            joinpath(CONFIG["aca_data"]["ACA_DATA_DIR"], aca_subdir, "$(layer_name).geojson")
        )

        @debug "Reprojecting ACA to target CRS"
        @time target_polys = GDF.reproject(target_polys, GI.crs(target_polys), EPSG_7844)

        @debug "Building STR tree for ACA"
        tree = STRT.STRtree(target_polys.geometry)
        reg_poly_idx = unique(vcat(STRT.query.(Ref(tree), management_zones.SHAPE)...))
        target_polys = target_polys[reg_poly_idx, :]

        # Standardize text
        target_polys.class .= lowercase.(
            replace.(target_polys.class, " " => "_", "/" => "_")
        )

        target_polys.class_id = map(
            x -> Symbol(x) in keys(id_lookup) ? getindex(id_lookup, Symbol(x)) : 0,
            target_polys.class
        )

        @debug "Caching reprojected ACA $(layer_name) data"
        GDF.write(aca_cache_fn, target_polys; crs=EPSG_7844)
    end

    # Rebuild query tree from filtered/standardized polygons for per-region queries
    tree = STRT.STRtree(target_polys.geometry)

    @info "Prepping hybrid $(layer_name) data"
    for reg in REGIONS
        fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_hybrid_$(layer_name).tif")
        if isfile(fn)
            @info "Skipping $reg as file already exists..."
            continue
        end

        force_gc_cleanup()

        reg_idx = region_geom_index(management_zones, reg; col=:AREA_DESCR)
        r = management_zones[reg_idx, :]

        @info "Cropping raster to management zone $reg"
        @time cropped_gbr10 = Rasters.trim(
            Rasters.mask(
                Rasters.crop(mpa_data; to=r.SHAPE);
                with=r.SHAPE
            )
        )

        # Select classes of interest and reproject to target CRS to ensure alignment
        @info "Ensuring reprojection is EPSG:7844"
        tmp_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_$(layer_name)_tmp.tif")
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
        # Use an anonymous lambda rather than `!=(0)` (a `Base.Fix` closure): ConstructionBase
        # v1.6.0 (latest release) cannot reconstruct `Base.Fix` on Julia 1.12+, which breaks
        # writing this lazy broadcast to disk. Anonymous lambdas take a different, working
        # code path in ConstructionBase/Flatten. See:
        # https://github.com/JuliaObjects/ConstructionBase.jl
        @time write_cog(fn, Int8.(ifelse.((x -> x != 0).(cropped_gbr10), cropped_gbr10, cropped_aca)))
        rm(tmp_fn)
    end

    return nothing
end
