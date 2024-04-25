using Rasters, YAXArrays, Rasters

import GeoDataFrames as GDF
import GeoFormatTypes as GFT
import ArchGDAL as AG

using GLMakie, GeoMakie


# tsv_data = joinpath(MPA_DATA_DIR, "Townsville-Whitsunday/SDB_AUS_TownsvilleWhitsunday_EOMAP_20180603_20180903_10m_MSL_geotiff.tif")
cairns_data = joinpath(MPA_DATA_DIR, "Cairns-Cooktown/SDB_AUS_CairnsCooktown_EOMAP_20180514_20181016_10m_MSL_geotiff.tif")
cairns_slope_data = joinpath(MPA_DATA_DIR, "slope/Cairns-Cooktown/slope_cnscook.tif")

moore_dpkg = "C:/Users/tiwanaga/development/ADRIA_data/Moore_2023-09-14_rc/site_data/Moore_2023-09-14_rc.gpkg"
gbr_dpkg = "C:/Users/tiwanaga/development/RME/rme_ml_2023_03_30b/data_files/region/reefmod_gbr.gpkg"

cairns_bathy = read(Raster(cairns_data, mappedcrs=EPSG(4326)))
cairns_slope = read(Raster(cairns_slope_data, mappedcrs=EPSG(4326)))
moore_gdf = GDF.read(moore_dpkg)

# crs(cairns_bathy)

# Reproject to match raster data (we're just assuming WGS84/EPSG:4326 is correct!)
moore_gdf.geom = GDF.reproject(moore_gdf.geom, GFT.EPSG(4326), GFT.EPSG(32755), order=:trad)

moore_bathy = trim(mask(cairns_bathy; with=moore_gdf.geom); pad=10)

f = plot(moore_bathy)
f.axis.title = "Moore Reef Cluster (Bathy)"
f.axis.xlabel = "Eastings"
f.axis.ylabel = "Northings"

moore_slope = trim(mask(cairns_slope; with=moore_gdf.geom); pad=10)
f = plot(moore_slope)
f.axis.title = "Moore Reef Cluster (Slope)"
f.axis.xlabel = "Eastings"
f.axis.ylabel = "Northings"


gbr_gdf = GDF.read(gbr_dpkg)
gbr_gdf.geom = GDF.reproject(gbr_gdf.geom, GFT.EPSG(4326), GFT.EPSG(32755), order=:trad)

cairns_region = trim(mask(cairns_bathy; with=gbr_gdf.geom); pad=100)

# Resampling should just be for visualization purposes!
cairns_low = resample(cairns_region, 100)
f = plot(cairns_low)
f.axis.title = "Cairns Region (Bathy)"
f.axis.xlabel = "Eastings"
f.axis.ylabel = "Northings"


cairns_regional_slope = trim(mask(cairns_slope; with=gbr_gdf.geom); pad=100)

cairns_slope_low = resample(cairns_regional_slope, 50)  # resample to lower resolution
f = plot(cairns_slope_low)
f.axis.title = "Cairns Region (Slope)"
f.axis.xlabel = "Eastings"
f.axis.ylabel = "Northings"
