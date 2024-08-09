using Statistics, StatsBase
using Distributed
using TOML
using GLMakie, GeoMakie
using ProgressMeter


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
            using TOML
            using Statistics, StatsBase
            using Glob
            using ProgressMeter
            using FLoops

            using Rasters
            using NCDatasets
            using GeoInterface
            using Distances

            using DataFrames
            import GeoDataFrames as GDF
            import GeoFormatTypes as GFT
            import ArchGDAL as AG

            using ImageCore: Gray
            using ImageFiltering
            using ImageContrastAdjustment: adjust_histogram, LinearStretching
            using ImageMorphology: label_components, component_centroids
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
        global GDA2020_crs = crs(regions_GDA2020[1, :geometry])

        global EPSG_4326 = EPSG(4326)  # Web mercator
        global EPSG_7844 = EPSG(7844)  # GDA2020 in degree projection
        global EPSG_7856 = EPSG(7856)  # GDA2020 in meter projection

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
        global WAVE_REGION_CRS = Dict(
            "Townsville-Whitsunday" => "+proj=utm +zone=55 +south +datum=WGS84",
            "Cairns-Cooktown" => "+proj=utm +zone=55 +south +datum=WGS84",
            "Mackay-Capricorn" => "+proj=utm +zone=56 +south +datum=WGS84",
            "FarNorthern" => "+proj=utm +zone=54 +south +datum=WGS84"
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
        ygridwidth=0.5,
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
function force_gc_cleanup()::Nothing
    sleep(1)  # Wait a little bit to ensure garbage sweep has occurred
    GC.gc()

    return nothing
end
