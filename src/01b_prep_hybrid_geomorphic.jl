"""
Creating a hybrid data product combining UQ-GBRMPA dataset with ACA to gap fill missing
(in-shore) areas of the GBR.

This script creates the hybridization of the two data products for geomorphic data.
"""

include("common.jl")

management_zones = GDF.read(
    joinpath(CONFIG["gda2020_data"]["GDA2020_DATA_DIR"], "Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg")
)

mpa_geomorphic_data = Raster(
    joinpath(CONFIG["mpa_data"]["MPA_DATA_DIR"], "geomorphic", "GBR10 GBRMP Geomorphic.tif");
    lazy=true
)

target_polys = GDF.read(
    joinpath(CONFIG["aca_data"]["ACA_DATA_DIR"], "Geomorphic-Map", "geomorphic.geojson")
)

# Reproject ACA data to target CRS
target_polys = GDF.reproject(target_polys, GI.crs(target_polys), EPSG_7844)

tree = STRT.STRtree(target_polys.geometry)
reg_poly_idx = unique(vcat(STRT.query.(Ref(tree), management_zones.SHAPE)...))
target_polys = target_polys[reg_poly_idx, :]

# Standardize text
target_polys.class .= lowercase.(replace.(target_polys.class, " " => "_"))
target_polys.class_id = map(x -> Symbol(x) in keys(MPA_GEOMORPHIC_IDS) ? getindex(MPA_GEOMORPHIC_IDS, Symbol(x)) : 0, target_polys.class)

@info "Prepping hybrid geomorphic data"
for reg in REGIONS
    fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_hybrid_geomorphic.tif")
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
        cropped_mz;
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

GC.gc()
