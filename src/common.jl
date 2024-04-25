using Distributed
using TOML

try
    global CONFIG = TOML.parsefile(".config.toml")
catch err
    if occursin("No such file", string(err))
        @info "Configuration file `.config.toml` not found. See README for instructions."
    end

    rethrow(err)
end

if nworkers() == 1 && (CONFIG["processing"]["N_PROCS"] > 1)
    addprocs(CONFIG["processing"]["N_PROCS"]; dir=@__DIR__)
end

@everywhere begin
    @eval begin
        using TOML
        using Statistics, StatsBase
        using Glob
        using ProgressMeter
        using FLoops

        using Rasters
        using NCDatasets

        import GeoDataFrames as GDF
        import GeoFormatTypes as GFT
        import ArchGDAL as AG

        using ImageCore: Gray
        using ImageFiltering
        using ImageContrastAdjustment: adjust_histogram, LinearStretching
        using ImageMorphology: label_components, component_centroids
    end

    global FIG_DIR = "../figs/"
    global OUTPUT_DIR = "../outputs/"
    global QGIS_DIR = "../qgis/"

    global CONFIG = TOML.parsefile(".config.toml")
    global MPA_DATA_DIR = CONFIG["mpa_data"]["MPA_DATA_DIR"]
    global ALLEN_ATLAS_DIR = CONFIG["aca_data"]["ALLEN_ATLAS_DIR"]

    # Folder names (TODO: Generalize)
    global REGIONS = String[
        "Townsville-Whitsunday",
        "Cairns-Cooktown",
        "Mackay-Capricorn",
        "FarNorthern",
    ]

    # Manually extracted from Raster Attribute Table(s)
    global FLAT_IDS = [13, 14, 23]  # Inner Reef Flat, Outer Reef Flat, Plateau
    global SLOPE_IDS = [21, 22, 24]  # Sheltered Reef Slope, Reef Slope, Back Reef Slope
    global BENTHIC_IDS = [0x0d, 0x0f]  # Rock, Coral/Algae

    # Known Proj4 strings for each GBRMPA zone
    global WAVE_REGION_CRS = Dict(
        "Townsville-Whitsunday" => "+proj=utm +zone=55 +south +datum=WGS84",
        "Cairns-Cooktown" => "+proj=utm +zone=55 +south +datum=WGS84",
        "Mackay-Capricorn" => "+proj=utm +zone=56 +south +datum=WGS84",
        "FarNorthern"=>"+proj=utm +zone=54 +south +datum=WGS84"
    )
end
