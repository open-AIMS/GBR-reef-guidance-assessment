"""
Creating a hybrid data product combining UQ-GBRMPA dataset with ACA to gap fill missing
(in-shore) areas of the GBR.

This script creates the hybridization of the two data products for benthic data.
"""

include("common.jl")

if !@isdefined(management_zones)
    management_zones = GDF.read(
        joinpath(CONFIG["gda2020_data"]["GDA2020_DATA_DIR"], "Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg")
    )
end

mpa_benthic_data = Raster(
    joinpath(CONFIG["mpa_data"]["MPA_DATA_DIR"], "benthic", "GBR10 GBRMP Benthic.tif");
    lazy=true,
    missingval=0
)

aca_benthic_data = GDF.read(
    joinpath(CONFIG["aca_data"]["ACA_DATA_DIR"], "Benthic-Map", "benthic.geojson")
)

target_polys = aca_benthic_data[aca_benthic_data.class.∈Ref(ACA_BENTHIC_IDS), :]
aca_benthic_data = nothing

# Reproject ACA data to target CRS
target_polys = GDF.reproject(target_polys, GI.crs(target_polys), EPSG_7844)

tree = STRT.STRtree(target_polys.geometry)
reg_poly_idx = vcat(STRT.query.(Ref(tree), management_zones.SHAPE)...)
target_polys = target_polys[reg_poly_idx, :]

# Standardize text
target_polys.class .= lowercase.(replace.(target_polys.class, " " => "_", "/" => "_"))
target_polys.class_id = map(x -> Symbol(x) in keys(MPA_BENTHIC_IDS) ? getindex(MPA_BENTHIC_IDS, Symbol(x)) : 0, target_polys.class)

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
    @time cropped_mz = Rasters.trim(
        Rasters.mask(
            Rasters.crop(mpa_benthic_data; to=r.SHAPE);
            with=r.SHAPE
        )
    )

    @info "Marking valid areas"
    cropped_mz = read(cropped_mz .* Bool.(cropped_mz .∈ Ref(values(MPA_BENTHIC_IDS))))

    # Select benthic classes of interest and reproject to target CRS to ensure alignment
    @info "Ensuring reprojection is EPSG:7844"
    tmp_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_benthic_tmp.tif")
    @time Rasters.resample(
        cropped_mz;
        crs=EPSG_7844,
        filename=tmp_fn
    )
    cropped_mz = Raster(tmp_fn; lazy=true, missingval=0)

    tree = STRT.STRtree(target_polys.geometry)
    reg_poly_idx = vcat(STRT.query.(Ref(tree), r.SHAPE)...)
    reg_polys = target_polys[reg_poly_idx, :]

    @info "Rasterizing polygons"
    @time d = Rasters.rasterize(
        maximum,
        reg_polys;
        to=cropped_mz,
        fill=:class_id,
        missingval=0
    )

    # target_benthic = crop_to_region(
    #     raw_benthic_fn,
    #     regions_4326[reg_idx_4326, :geometry],
    #     criteria_paths[:Benthic]
    # )
    # resample_and_write(
    #     target_benthic,
    #     bathy_gda2020,
    #     criteria_paths[:Benthic];
    #     method=:near
    # )
    # target_benthic = nothing
    # force_gc_cleanup()

    @info "Writing hybrid data for $reg"
    write_cog(fn, Raster(cropped_mz; data=sparse(Int8.(cropped_mz .| d))))
    rm(tmp_fn)
end

mpa_benthic_data = nothing
target_polys = nothing

# Clear remnants of this processing
GC.gc()
