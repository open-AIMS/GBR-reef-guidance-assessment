using Distributed

const DATA_DIR = "D:/development/ADRIA_data/spatial_datasets/Bathy data 10m"
const RESULT_DIR = "../figs/"
const QGIS_DIR = "../qgis/"

# Folder names (TODO: Generalize)
REGIONS = String[
    "Townsville-Whitsunday",
    "Cairns-Cooktown",
    "Mackay-Capricorn",
    "FarNorthern",
]

# Manually extracted from Raster Attribute Table(s)
FLAT_IDS = [13, 14, 23]  # Inner Reef Flat, Outer Reef Flat, Plateau
SLOPE_IDS = [21, 22, 24]  # Sheltered Reef Slope, Reef Slope, Back Reef Slope
BENTHIC_IDS = [0x0d, 0x0f]  # Rock, Coral/Algae

if nworkers() == 1
    addprocs(2)
    # addprocs(4, exeflags = "--project=$(Base.active_project())")

    @everywhere begin
        const DATA_DIR = "D:/development/ADRIA_data/spatial_datasets/Bathy data 10m"
        const RESULT_DIR = "../figs/"
        const QGIS_DIR = "../qgis/"

        # Folder names (TODO: Generalize)
        REGIONS = String[
            "Townsville-Whitsunday",
            "Cairns-Cooktown",
            "Mackay-Capricorn",
            "FarNorthern",
        ]

        # Manually extracted from Raster Attribute Table(s)
        FLAT_IDS = [13, 14, 23]  # Inner Reef Flat, Outer Reef Flat, Plateau
        SLOPE_IDS = [21, 22, 24]  # Sheltered Reef Slope, Reef Slope, Back Reef Slope
        BENTHIC_IDS = [0x0d, 0x0f]  # Rock, Coral/Algae

        using Rasters

        import GeoDataFrames as GDF
        import GeoFormatTypes as GFT
        import ArchGDAL as AG

        using ImageCore: Gray
        using ImageFiltering
        using ImageContrastAdjustment: adjust_histogram, LinearStretching
        using ImageMorphology: label_components, component_centroids

        using Statistics, StatsBase
        using Glob
        using ProgressMeter
    end
end