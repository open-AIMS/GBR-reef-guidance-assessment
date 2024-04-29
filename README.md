# Bathy10m

Analyses to support workshop discussions.

## Project Layout

Assumes `src` is the project root. Each file in `src` is expected to be run in order.

```code
GBR-FeatureAnalysis/
├─ src/           # Analysis
├─ outputs/       # Intermediate data files
├─ figs/          # Figures
├─ .gitignore
├─ Project.toml   # Julia project spec
├─ LICENSE.md
├─ README.md      # this file
```

In addition to the above, spatial data should be stored outside the repository and its
location defined in a `.config.toml` file placed within `src`.

```TOML
[processing]
N_PROCS = 2  # Number of cores to use for multi-processing steps

[mpa_data]
MPA_DATA_DIR = "path to GBR data"  # location of GBR datasets

[aca_data]
ALLEN_ATLAS_DIR = "path to Allen Atlas data"  # location of Allen Atlas datasets
```

### Data layout

Expected data directory layout:

Sub-directory names should be consistent and match.
Benthic habitat raster (inside `benthic`) covers whole of GBR so no sub-directories are necessary.
Similarly, `geomorphic` holds the whole-of-GBR geomorphic zonation raster

`zones` holds GBRMPA zone layers in geojson format.

`features` holds the GBRMPA GBR-wide feature set.

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
├───waves
│   ├───Cairns-Cooktown
│   ├───FarNorthern
│   ├───Mackay-Capricorn
│   └───Townsville-Whitsunday
└───zones
```

## Scripts

All scripts assume `src` is the root. References to external directories may be specified
as relative to the `src` directory.

```bash
$ cd src
$ julia --project=..
```

Scripts are labelled by their expected run order, and are written to be as stand-alone as
possible. It should be possible to run one script, so long as other scripts earlier in the
indicated order have been run previously.

e.g., Script 3 could be run after script 1 and 2, so long as 1 and 2 were run at some point
previously.

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

#### Projections

Wave data provided can span across two UTM zones.
Saving processed data as geotiffs will fail due to bounds checking; the GDAL implementation
will refuse to write to disk if the data falls outside of the indicated projection bounds.

This may be why the wave data (provided in netCDF format) does not include a CRS.

For example, most of the Townsville-Whitsunday region falls under UTM Zone 55S extents
(EPSG: 32755).

Latitude: 1116915.04, 10000000.0
Longitude: 166021.44, 833978.56

Units are meters in Northings/Eastings. Simply, the greater the Longitude (Eastings) value,
the further east the extent. The greater the Latitude (Northings) value, the further north
the extent.

The bathymetry data for the same region spans across:

Latitude: 7709845.2670042375, 8045915.2670042375
Longitude: 394844.2269136696, 875934.2269136696

The wave data spans across:

Latitude: 7709850.2670042375, 8045910.2670042375
Longitude: 394849.2269136696, 875929.2269136696

After trimming to bathymetry extents:

Latitude: 7729280.2670042375, 8040970.2670042375
Longitude: 401619.2269136696, 859799.2269136696

Note the east-most longitude is outside the nominal bounds for UTM Zone 55S.

As a tentative workaround, it is assumed that the embedded extents in the wave netCDFs are
incorrect. Extents of the wave data typically match the size/shape of bathymetry data
(i.e., the number of rows/columns are the same). When missing data is cropped away, the
remaining wave data is well within the bathymetry bounds. Where the sizes do not match,
the bounds of the wave data is extended so it does (the additional rows/columns are filled
with values denoting `missing` data).

The bathymetry data structure is then copied (so metadata on its extent, projection, etc.
are retained), and finally, the wave data is copied across.

## Presentation

Recommended Colors

- Flats: #7570b3  (purple)
- Slopes: #1b9e77  (green)

https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
