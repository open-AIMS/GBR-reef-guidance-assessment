# Bathy10m

Analyses to support workshop discussions.

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

Once the layers have been added, create an MBTiles file with:

Processing Toolbox -> Raster Tools -> Generate XYZ tiles (MBTiles)

## Data Sources

Benthic Habitat Layer:
- Great Barrier Reef 10m Grid (GBR10) GBRMP Benthic
  Great Barrier Reef Marine Park Authority
  https://gbrmpa.maps.arcgis.com/home/item.html?id=d1c58d71667d490ba650c8fd07d6f7ee
  https://metadata.imas.utas.edu.au/geonetwork/srv/eng/catalog.search#/metadata/492a87d95e8243728486718e7aed02a8


Bathymetry 10m Grid:
- https://gbrmpa.maps.arcgis.com/home/item.html?id=f644f02ec646496eb5d31ad4f9d0fc64

Slope 10m Grid:

To be confirmed.

GBRMPA Features:
- https://data.gov.au/dataset/ds-dga-51199513-98fa-46e6-b766-8e1e1c896869/details

Geomorphic:
- https://gbrmpa.maps.arcgis.com/home/item.html?id=93fd689452e44e74801845b7935c54c4

GBRMPA Zones:
- https://geoportal.gbrmpa.gov.au/datasets/GBRMPA::great-barrier-reef-marine-park-zoning/explore
- https://geoportal.gbrmpa.gov.au/datasets/GBRMPA::management-areas-of-the-great-barrier-reef-marine-park/explore?location=-17.583829%2C150.586624%2C6.38


GBRMPA Bioregions:
- https://geoportal.gbrmpa.gov.au/datasets/GBRMPA::reef-marine-bioregions-of-the-great-barrier-reef/about

Note: Bathymetry and slope data obtained via M. Poutinen.

## Data layout

Expected data directory layout:

Sub-directory names should be consistent and match.
Benthic habitat raster (inside `benthic`) covers whole of GBR so no sub-directories are necessary.
Similarly, `geomorphic` holds the whole-of-GBR geomorphic zonation raster

`zones` holds GBRMPA zone layers in geojson format.

`features` holds the GBRMPA GBR-wide feature set.

```bash
DATA_DIR
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
```

Recommended Colors

- Flats: #7570b3  (purple)
- Slopes: #1b9e77  (green)

https://colorbrewer2.org/#type=qualitative&scheme=Dark2&n=3
