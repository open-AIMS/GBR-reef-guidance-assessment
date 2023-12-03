# Bathy10m

Analyses to support workshop discussions.

## Scripts

All scripts assume `src` is the root. References to external directories may be specified
as relative to the `src` directory.

```bash
$ cd src
$ julia --project=..
```

Scripts are labelled by their expected run order.


## Data Sources

Benthic Habitat Layer:
    Great Barrier Reef 10m Grid (GBR10) GBRMP Benthic
    Great Barrier Reef Marine Park Authority
    https://gbrmpa.maps.arcgis.com/home/item.html?id=d1c58d71667d490ba650c8fd07d6f7ee
    https://metadata.imas.utas.edu.au/geonetwork/srv/eng/catalog.search#/metadata/492a87d95e8243728486718e7aed02a8


Bathymetry 10m Grid:
    https://gbrmpa.maps.arcgis.com/home/item.html?id=f644f02ec646496eb5d31ad4f9d0fc64

Slope 10m Grid:
    To be confirmed.

Geomorphic:
    Currently unused.

    https://gbrmpa.maps.arcgis.com/home/item.html?id=93fd689452e44e74801845b7935c54c4


Bathymetry and slope data obtained via M. Poutinen.

## Data layout

Expected data directory layout:

Sub-directory names should be consistent and match.
Benthic habitat raster covers whole of GBR so no sub-directories are necessary.

```bash
DATA_DIR
├───bathy
│   ├───Cairns-Cooktown
│   ├───FarNorthern
│   ├───Mackay-Capricorn
│   └───Townsville-Whitsunday
├───benthic
└───slope
    ├───Cairns-Cooktown
    ├───FarNorthern
    ├───Mackay-Capricorn
    └───Townsville-Whitsunday
```