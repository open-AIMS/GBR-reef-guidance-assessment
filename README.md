# GBR-reef-guidance-assessment

Analyses to support workshop discussions for pilot deployment program.

## Project Layout

Assumes `src` is the project root. Each file in `src` is expected to be run in order.

```code
GBR-reef-guidance-assessment/
├─ src/           # Analysis
├─ outputs/       # Intermediate data files
├─ figs/          # Figures
├─ .gitignore
├─ Project.toml   # Julia project spec
├─ LICENSE.md
└─ README.md      # this file
```

In addition to the above, spatial data should be stored outside the repository and its
location defined in a `.config.toml` file placed within `src`.

```TOML
[processing]
N_PROCS = 2  # Number of cores to use for multi-processing steps

[mpa_data]
MPA_DATA_DIR = "path to GBR data"  # location of GBR datasets

[aca_data]
ACA_DATA_DIR = "path to Allen Atlas data"  # location of Allen Atlas datasets

[wave_data]
WAVE_DATA_DIR = "path to wave data"  # location of wave NetCDF datasets

[gda2020_data]
GDA2020_DATA_DIR = "path to GDA2020 vector data"  # location of GBRMPA region feature datasets in EPSG:7844 (GDA2020)

[rugosity_data]
RUG_DATA_DIR = "path to Rugosity raster data"  # location of Rugosity data

[ports_data]
PORT_DATA_DIR = "path to QLD Ports vector data"  # Location of QLD Ports data.
```

### Data layout

Expected data directory layout:

Separate directories are used for GBRMPA (MPA), Allen Coral Atlas (ACA), wave, GDA2020 and Rugosity datasets.
Sub-directory names should be consistent and match.
Data directories can currently be accessed via `AIMS-Decision Support Sharepoint/Documents/General/GBR Spatial Datasets`.

MPA_DATA_DIR : contains raster data at whole-GBR and GBRMPA-management-region scales.
- `zones` holds GBRMPA zone layers in geojson format.
- `features` holds the GBRMPA GBR-wide feature set.
- Sharepoint folder name : `GBR-Bathy10m`.

ACA_DATA_DIR : contains raster and vector data at whole-GBR scale.
- Sharepoint folder name : `AllenAtlas_GBR-20231118074407`.

WAVE_DATA_DIR : contains wave data in NetCDF format at the scale of GBRMPA-management-regions.
- Sharepoint folder name : `Wave-Data-for-PDP`.

GDA2020_DATA_DIR : contains GBRMPA management and zoning geopackage data in EPSG:7844 (GDA2020).
- Sharepoint folder name: `GDA2020-Data-for-PDP`.

RUG_DATA_DIR : contains Rugosity raster data provided by Ben Radford for Townsville-Whitsunday region.
- Sharepoint folder name : `GBR-Rugosity_Radford`.

PORT_DATA_DIR : contains QLD Port location vector data provided by Marji Poutinen.
- Sharepoint folder name : `QLD_ports_mercator_via_MP`.

```bash
MPA_DATA_DIR
├───bathy
│   ├───Cairns-Cooktown
│   ├───FarNorthern
│   ├───Mackay-Capricorn
│   └───Townsville-Whitsunday
├───benthic
├───features
├───geomorphic
├───slope
│   ├───Cairns-Cooktown
│   ├───FarNorthern
│   ├───Mackay-Capricorn
│   └───Townsville-Whitsunday
└───zones

ACA_DATA_DIR
├───Bathymetry---composite-depth
│   └───Raster data
├───Benthic-Map
│   └───Vector data
├───boundary
│   └───Vector data
├───Geomorphic-Map
│   └───Vector data
├───Reef-Extent
│   └───Vector data
└───Turbidity-Q3-2023
    └───Raster data

GDA2020_DATA_DIR
├───management_region_features
└───marine_park_zoning

WAVE_DATA_DIR
├───Hs
│   ├───Cairns-Cooktown
│   ├───FarNorthern
│   ├───Mackay-Capricorn
│   └───Townsville-Whitsunday
└───Tp
    ├───Cairns-Cooktown
    ├───FarNorthern
    ├───Mackay-Capricorn
    └───Townsville-Whitsunday

RUG_DATA_DIR
└───Townsville-Whitsunday_Rugosity Raster

PORT_DATA_DIR
└───QLD_Ports Shapefile
```

## Scripts

All scripts assume `src` is the root. References to external directories may be specified
as relative to the `src` directory.

```bash
$ cd src
$ julia --project=..
```

Scripts are labelled by their expected run order for GBRMPA (MPA) and Allen Coral Atlas (ACA),
and are written to be as stand-alone as possible.
It should be possible to run one script, so long as other scripts earlier in the
indicated order have been run previously.
e.g., Script 3 could be run after script 1 and 2, so long as 1 and 2 were run at some point
previously.
(Note: Running `1b_*.jl` requires the user to have raw MPA bathymetry input data,
as it is required for processing of wave data.)
(Note: Running `1a_*.jl` requires the user to have raw ACA turbidity input data,
as it is also required for MPA analyses.)

- `1_*.jl` : Separate data into regions to reduce computational requirements. Ensure all data
used in later steps are in EPSG:7844/GDA2020.
- `2_*.jl` : Filter raster data into cells that meet selected criteria and calculate the
proportion of suitability in the hectare surrounding each cell.
- `3_*.jl` : Count the number of cells that have a surrounding suitability >= 0.95.

## Manual steps

In QGIS, once the raster layers have been added, create an MBTiles file with:

Processing Toolbox -> Raster Tools -> Generate XYZ tiles (MBTiles)

## Data Sources

### Benthic Habitat Layer

Great Barrier Reef 10m Grid (GBR10) GBRMP Benthic
Great Barrier Reef Marine Park Authority
https://gbrmpa.maps.arcgis.com/home/item.html?id=d1c58d71667d490ba650c8fd07d6f7ee
https://metadata.imas.utas.edu.au/geonetwork/srv/eng/catalog.search#/metadata/492a87d95e8243728486718e7aed02a8

### Bathymetry 10m Grid

https://gbrmpa.maps.arcgis.com/home/item.html?id=f644f02ec646496eb5d31ad4f9d0fc64

### Slope 10m Grid

Calculated by Dr M. Puotinen based on the bathymetry.

### GBRMPA Features

https://data.gov.au/dataset/ds-dga-51199513-98fa-46e6-b766-8e1e1c896869/details

### Geomorphic

https://gbrmpa.maps.arcgis.com/home/item.html?id=93fd689452e44e74801845b7935c54c4

### GBRMPA Zones

https://geoportal.gbrmpa.gov.au/datasets/GBRMPA::management-areas-of-the-great-barrier-reef-marine-park/explore?location=-17.583829%2C150.586624%2C6.38
https://geoportal.gbrmpa.gov.au/datasets/GBRMPA::great-barrier-reef-marine-park-zoning/explore

### GBRMPA Bioregions

https://geoportal.gbrmpa.gov.au/datasets/GBRMPA::reef-marine-bioregions-of-the-great-barrier-reef/about

### Wave data

Callaghan, David (2023). Great Barrier Reef non-cyclonic and on-reef wave model predictions.
The University of Queensland.
Data Collection.
https://doi.org/10.48610/8246441
https://espace.library.uq.edu.au/view/UQ:8246441

#### Notes

`Hs` (Significant wave height) : average height of the top 1/3 highest waves in a section
of ocean; in metres Gives a general idea of sea state, which is the general roughness
(or not) of a part of the sea A reasonable proxy for the potential for wave impacts on
structures like reefs.
Often this might be the only data available.

`Tp` (Peak wave period) : Peak wave energy period, in seconds.
A low `Tp` means a locally generated wind wave.
A high `Tp` means swell generated far away that has propagated to the location of interest.
If you need to know whether the sea state is generated by local winds or from a distant storm.
`Tp` < 6 s : Mixed or wind wave environments
`Tp` > 8 s : Swell wave environments

`Tp` is used to constrain locations by suitable deployment conditions.
Deployment conditions are dependent on the type of vessel. Small tenders are assumed to be
used for pilot deployment scenarios for 2025.

`Hs` and `Tp` data from the 90th percentile of values are used (`Hs90` and `Tp90`).

### ACA data

https://www.allencoralatlas.org/

### Rugosity data

Provided by Ben Radford for Townsville-Whitsunday region.

## Resolution

`MPA raster data` (Bathymetry, Benthic, Geomorphic, Slope) : 10 x 10m pixel size

`ACA raster data` (Bathymetry and Turbidity) : 10 x 10m pixel size

`Waves NetCDF data` (Hs and Tp) : 10 x 10m pixel size

### Projections

After `1*_.jl` (pre-processing and reprojecting data) is complete all raster and vector data
in `*_OUTPUT_DIR`, and outputs from `2*_.jl` and `3*_.jl` will be in EPSG:7844 GDA2020.

#### MPA Data

Wave data provided can span across two UTM zones.
Saving processed data as geotiffs will fail due to bounds checking; the GDAL implementation
will refuse to write to disk if the data falls outside of the indicated projection bounds.

This may be why the wave data (provided in netCDF format) does not include a CRS.

For example, most of the Townsville-Whitsunday region falls under UTM Zone 55S extents
(EPSG:32755).

- Latitude: 1116915.04, 10000000.0
- Longitude: 166021.44, 833978.56

Units are meters in Northings/Eastings. Simply, the greater the Longitude (Eastings) value,
the further east the extent. The greater the Latitude (Northings) value, the further north
the extent.

The bathymetry data for the same region spans across:

- Latitude: 7709845.2670042375, 8045915.2670042375
- Longitude: 394844.2269136696, 875934.2269136696

The wave data spans across:

- Latitude: 7709850.2670042375, 8045910.2670042375
- Longitude: 394849.2269136696, 875929.2269136696

After trimming to bathymetry extents:

- Latitude: 7729280.2670042375, 8040970.2670042375
- Longitude: 401619.2269136696, 859799.2269136696

**Note the east-most longitude is outside the nominal bounds for UTM Zone 55S.**

As a tentative workaround, it is assumed that the embedded extents in the wave netCDFs are
incorrect. Extents of the wave data typically match the size/shape of bathymetry data
(i.e., the number of rows/columns are the same). When missing data is cropped away, the
remaining wave data is well within the bathymetry bounds. Where the sizes do not match,
the bounds of the wave data is extended so it does (the additional rows/columns are filled
with values denoting `missing` data).

The bathymetry data structure is then copied (so metadata on its extent, projection, etc.
are retained), and finally, the wave data is copied across.

All data in `1a_prep_MPA.jl` are projected to EPSG:7844/GDA2020 prior to further analysis.

#### ACA Data

ACA raster data is in EPSG:4326/WGS84.
Wave data preparation in `1b_prep_aca.jl` requires
MPA bathymetric data to be available due to inconsistent reference systems/units
(see 'Projections - MPA Data' section for more details).

All data in `1b_prep_aca.jl` are projected to EPSG:7844/GDA2020 prior to further analysis.

## Presentation

Recommended Colors

- Flats: #7570b3  (purple)
- Slopes: #1b9e77  (green)

https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
