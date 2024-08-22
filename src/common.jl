using Statistics, StatsBase
using Distributed
using TOML
using GLMakie, GeoMakie
using ProgressMeter
import GeometryOps as GO


try
    global CONFIG = TOML.parsefile(".config.toml")
catch err
    if occursin("No such file", string(err))
        @info "Configuration file `.config.toml` not found. See README for instructions."
    end

    rethrow(err)
end

if nworkers() < CONFIG["processing"]["N_PROCS"]
    addprocs(CONFIG["processing"]["N_PROCS"]; dir=@__DIR__)

    @everywhere begin
        @eval begin
            using Statistics, StatsBase
            using ProgressMeter
            using TOML
            using Glob

            using
                FLoops,
                ThreadsX

            using
                Rasters,
                NCDatasets,
                LibGEOS

            import GeoInterface as GI
            import GeometryOps as GO
            import GeoFormatTypes as GFT
            import ArchGDAL as AG

            using Distances

            using DataFrames
            import GeoDataFrames as GDF
            import GeoParquet as GP

            using Images, ImageFiltering
            using ImageMorphology: label_components
        end

        FIG_DIR = "../figs/"
        global MPA_FIG_DIR = joinpath(FIG_DIR, "MPA")
        global ACA_FIG_DIR = joinpath(FIG_DIR, "ACA")

        QGIS_DIR = "../qgis/"
        global MPA_QGIS_DIR = joinpath(QGIS_DIR, "MPA")
        global ACA_QGIS_DIR = joinpath(QGIS_DIR, "ACA")

        OUTPUT_DIR = "../outputs/"
        global MPA_OUTPUT_DIR = joinpath(OUTPUT_DIR, "MPA")
        global ACA_OUTPUT_DIR = joinpath(OUTPUT_DIR, "ACA")

        global MPA_ANALYSIS_RESULTS = joinpath(MPA_OUTPUT_DIR, "analysis_results")
        global ACA_ANALYSIS_RESULTS = joinpath(ACA_OUTPUT_DIR, "analysis_results")

        CONFIG = TOML.parsefile(".config.toml")
        global MPA_DATA_DIR = CONFIG["mpa_data"]["MPA_DATA_DIR"]
        global ACA_DATA_DIR = CONFIG["aca_data"]["ACA_DATA_DIR"]
        global WAVE_DATA_DIR = CONFIG["wave_data"]["WAVE_DATA_DIR"]
        global GDA2020_DATA_DIR = CONFIG["gda2020_data"]["GDA2020_DATA_DIR"]
        global RUG_DATA_DIR = CONFIG["rugosity_data"]["RUG_DATA_DIR"]
        global PORT_DATA_DIR = CONFIG["ports_data"]["PORT_DATA_DIR"]

        regions_GDA2020_path = joinpath(
            GDA2020_DATA_DIR,
            "Great_Barrier_Reef_Marine_Park_Management_Areas_20_1685154518472315942.gpkg"
        )
        regions_GDA2020 = GDF.read(regions_GDA2020_path)
        rename!(regions_GDA2020, Dict(:SHAPE => :geometry))

        global EPSG_4326 = GFT.EPSG(4326)  # Web mercator
        global EPSG_7844 = GFT.EPSG(7844)  # GDA2020 in degrees
        global EPSG_9473 = GFT.EPSG(9473)  # GDA2020 in meter projection

        # Get polygon of management areas
        global REGION_PATH_4326 = joinpath(
            MPA_DATA_DIR,
            "zones",
            "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
        )

        # The reef features GDA94 dataset has added `Area_HA` column needed in 3*_.jl
        global REEF_PATH_GDA94 = joinpath(
            MPA_DATA_DIR,
            "features",
            "Great_Barrier_Reef_Features.shp"
        )

        # Folder names (TODO: Generalize)
        global REGIONS = String[
            "Townsville-Whitsunday",
            "Cairns-Cooktown",
            "Mackay-Capricorn",
            "FarNorthern"
        ]

        # GBRMPA IDs Manually extracted from Raster Attribute Table(s)
        global MPA_FLAT_IDS = [13, 14, 23]  # Inner Reef Flat, Outer Reef Flat, Plateau
        global MPA_SLOPE_IDS = [21, 22, 24]  # Sheltered Reef Slope, Reef Slope, Back Reef Slope
        global MPA_BENTHIC_IDS = [0x0d, 0x0f] # 0x0d = 13 = Rock # 0x0f = 15 = Coral/Algae

        # ACA IDs Manually extracted from Raster Attribute Table(s)
        global ACA_FLAT_IDS = [
            "Terrestrial Reef Flat", "Plateau", "Inner Reef Flat", "Outer Reef Flat"
        ]
        global ACA_SLOPE_IDS = ["Sheltered Reef Slope", "Back Reef Slope", "Reef Slope"]
        global ACA_BENTHIC_IDS = ["Coral/Algae", "Rock"]

        # Known Proj strings for each GBRMPA zone - may remove in later cleanup?
        global REGION_CRS_UTM = Dict(
            "Townsville-Whitsunday" => GFT.EPSG(32755),
            "Cairns-Cooktown" => GFT.EPSG(32755),
            "Mackay-Capricorn" => GFT.EPSG(32756),
            "FarNorthern" => GFT.EPSG(32754)
        )

        # GBRMPA zones to exclude from site selection
        global MPA_EXCLUSION_ZONES = ["Preservation Zone"]
    end
end

function plot_map(gdf::DataFrame; geom_col=:geometry, color=nothing)
    f = Figure(; size=(600, 900))
    ga = GeoAxis(
        f[1,1];
        dest="+proj=latlong +datum=WGS84",
        xlabel="Longitude",
        ylabel="Latitude",
        xticklabelpad=15,
        yticklabelpad=40,
        xticklabelsize=10,
        yticklabelsize=10,
        aspect=AxisAspect(0.75),
        xgridwidth=0.5,
        ygridwidth=0.5
    )

    plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))

    if !isnothing(color)
        poly!(ga, plottable, color=color)
    else
        poly!(ga, plottable)
    end

    # Need to auto-set limits explicitly, otherwise tick labels don't appear properly (?)
    # autolimits!(ga)
    xlims!(ga)
    ylims!(ga)

    display(f)

    return f
end
function plot_map!(ga::GeoAxis, gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing

    plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))

    if !isnothing(color)
        poly!(ga, plottable; color=color)
    else
        poly!(ga, plottable)
    end

    # Need to auto-set limits explicitly, otherwise tick labels don't appear properly (?)
    # autolimits!(ga)
    xlims!(ga)
    ylims!(ga)

    return nothing
end
function plot_map!(gdf::DataFrame; geom_col=:geometry, color=nothing)::Nothing
    ga = current_axis()
    plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))

    if !isnothing(color)
        poly!(ga, plottable; color=color)
    else
        poly!(ga, plottable)
    end

    # Need to auto-set limits explicitly, otherwise tick labels don't appear properly (?)
    # autolimits!(ga)
    xlims!(ga)
    ylims!(ga)

    return nothing
end

"""
    force_gc_cleanup()

Trigger garbage collection to free memory after clearing large datasets.
Not exactly best practice, but it works for very high memory workloads where data is
repeatedly loaded/unloaded.
"""
function force_gc_cleanup(; wait_time=1)::Nothing
    sleep(wait_time)  # Wait a little bit to ensure garbage sweep has occurred
    GC.gc()

    return nothing
end

"""
    port_buffer_mask(gdf::DataFrame, dist::Float64; unit::String="NM")

Create a masking buffer around indicated port locations.

# Arguments
- `gdf` : GeoDataFrame of port locations (given as long/lat points)
- `dist` : distance from port in degrees (deg), kilometers (km), or nautical miles (NM; default)
- `unit` : unit `dist` is in
"""
function port_buffer_mask(gdf::DataFrame, dist::Float64; unit::String="NM")
    # Determine conversion factor (nautical miles or kilometers)
    conv_factor = 1.0
    if unit == "NM"
        conv_factor = 60.0  # 60 NM = 1 degree
    elseif unit == "km"
        conv_factor = 111.0  # 111 km = 1 degree
    elseif unit != "deg"
        error("Unknown distance unit requested. Can only be one of `NM` or `km` or `deg`")
    end

    ports = gdf.geometry  # TODO: Replace with `GI.geometrycolumns()`

    # Make buffer around ports
    buffered_ports = GO.buffer.(ports, dist / conv_factor)

    # Combine all geoms into one
    port_mask = reduce((x1, x2) -> LibGEOS.union(x1, x2), buffered_ports)

    return port_mask
end
