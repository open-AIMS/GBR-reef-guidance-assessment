"""
Creating a hybrid data product combining UQ-GBRMPA dataset with ACA to gap fill missing
(in-shore) areas of the GBR.

This script creates the hybridization of the two data products for benthic data.
"""

include("common.jl")

management_zones = GDF.read(
    joinpath(CONFIG["gda2020_data"]["GDA2020_DATA_DIR"], "Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg")
)

mpa_benthic_data = Raster(
    joinpath(CONFIG["mpa_data"]["MPA_DATA_DIR"], "benthic", "GBR10 GBRMP Benthic.tif");
    lazy=true
)

aca_benthic_data = GDF.read(
    joinpath(CONFIG["aca_data"]["ACA_DATA_DIR"], "Benthic-Map", "benthic.geojson")
)

target_polys = aca_benthic_data[aca_benthic_data.class.∈Ref(ACA_BENTHIC_IDS), :]
aca_benthic_data = nothing

target_polys.class_id .= 0

for (aca, mpa) in zip(ACA_BENTHIC_IDS, MPA_BENTHIC_IDS)
    target_polys[target_polys.class.==aca, :class_id] .= mpa
end

# Reproject ACA data to target CRS
target_polys = GDF.reproject(target_polys, GI.crs(target_polys), EPSG_7844)

tree = STRT.STRtree(target_polys.geometry)
reg_poly_idx = vcat(STRT.query.(Ref(tree), management_zones.SHAPE)...)
target_polys = target_polys[reg_poly_idx, :]

@info "Prepping hybrid benthic data"
for reg in REGIONS
    fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_hybrid_benthic.tif")
    if isfile(fn)
        @info "Skipping $reg as file already exists..."
        continue
    end

    reg_idx = occursin.(reg[1:3], management_zones.AREA_DESCR)
    r = management_zones[reg_idx, :]

    @info "Cropping raster to management zone $reg"
    @time cropped_mz = Rasters.crop(mpa_benthic_data; to=r.SHAPE)

    # Select benthic classes of interest and reproject to target CRS to ensure alignment
    # if GI.crs(cropped_mz)
    @info "Ensuring reprojection is EPSG:7844"
    @time cropped_mz = Rasters.resample(
        Int8.(read(cropped_mz) .∈ Ref(MPA_BENTHIC_IDS));
        crs=EPSG_7844
    )

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

    @info "Writing hybrid data for $reg"
    Rasters.write(fn, Int8.(cropped_mz .| d))
end

# Clear remnants of this processing
GC.gc()
