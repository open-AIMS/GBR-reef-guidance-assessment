using Distributed

DATA_DIR = "C:/Users/rlippman/Documents/development/ADRIA_data/spatial_datasets/GBR-Bathy10m"
# "D:/development/ADRIA_data/spatial_datasets/Bathy data 10m"
RESULT_DIR = "../figs_Allen/"
QGIS_DIR = "../qgis/"

ALLEN_DATA_DIR = joinpath(DATA_DIR, "..", "AllenAtlas_GBR-20231118074407")
ALLEN_RESULT_DIR = joinpath(ALLEN_DATA_DIR, "outputs")
# LAYERS =  ["Bathymetry---composite-depth", "Turbidity-Q3-2023"] # these files to be converted to TIFs
#                 # "Benthic-Map", "Geomorphic-Map",              #readdir(DATA_DIR)    # Criteria = all files/folders in the directory

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
        DATA_DIR = "C:/Users/rlippman/Documents/development/ADRIA_data/spatial_datasets/GBR-Bathy10m"
        RESULT_DIR = "../figs/"
        QGIS_DIR = "../qgis/"

        ALLEN_DATA_DIR = joinpath(DATA_DIR, "..", "AllenAtlas_GBR-20231118074407")
        ALLEN_RESULT_DIR = joinpath(ALLEN_DATA_DIR, "outputs")
        # LAYERS =  ["Bathymetry---composite-depth", "Turbidity-Q3-2023"] # these files to be converted to TIFs
        #                 # "Benthic-Map", "Geomorphic-Map",              #readdir(DATA_DIR)    # Criteria = all files/folders in the directory

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

        BENTHIC_ALLEN = ["Rock", "Coral/Algae"]
        FLATS_ALLEN = ["Inner Reef Flat", "Outer Reef Flat", "Plateau"]           # identify flats
        SLOPES_ALLEN = ["Sheltered Reef Slope", "Reef Slope", "Back Reef Slope"]  # identify slopes
        ## Geomorphic zones unaccounted for
        # "Deep Lagoon" - UNWANTED
        # "Shallow Lagoon"
        # "Reef Crest" -  between the Reef Flat and the Reef Slope
        # "Terrestrial Reef Flat" - flat - wanted?
        # "Patch Reefs"


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