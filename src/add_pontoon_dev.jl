using ColorSchemes

include("common.jl")

boat_data = GDF.read(joinpath(MPA_DATA_DIR,"../operational_boating_facilities.gpkg"))

region_path = joinpath(
    MPA_DATA_DIR,
    "zones",
    "Management_Areas_of_the_Great_Barrier_Reef_Marine_Park.geojson"
)
region_features = GDF.read(region_path)

relevant_points = DataFrame(geometry = [], facility = [], buffer = [])

for point in eachrow(boat_data)
    if any(AG.intersects.([point.geom], region_features.geometry))
        push!(relevant_points, [point.geom, point.facility, missing])
    end
end

for point in eachrow(relevant_points)
    buffer = AG.buffer(point.geometry, 0.2)
    point.buffer = buffer
end

function _convert_plottable(gdf::Union{DataFrame,DataFrameRow}, geom_col::Symbol)
    local plottable
    try
        if gdf isa DataFrame
            plottable = GeoMakie.geo2basic(AG.forceto.(gdf[!, geom_col], AG.wkbPolygon))
        else
            plottable = GeoMakie.geo2basic(AG.forceto(gdf[geom_col], AG.wkbPolygon))
        end
    catch
        # Column is already in a plottable form, or some unrelated error occurred
        if gdf isa DataFrame
            plottable = gdf[:, geom_col]
        else
            plottable = [gdf[geom_col]]
        end
    end

    return plottable
end

"""
    plot_map(gdf::DataFrame; geom_col::Symbol=:geometry, color_by::Symbol)

Convenience plot function.

# Arguments
- `gdf` : GeoDataFrame
- `geom_col` : Column name holding geometries to plot
- `color_by` : Column name holding factor to color reefs by (e.g. :management_area)
- `buffer_col` : Column name holding geometries of buffer areas to plot
"""
function plot_map(gdf::DataFrame; geom_col::Symbol=:geometry, color_by::Symbol, buffer_col::Symbol)
    f = Figure(; size=(600, 900))
    ga = f[1, 1] = GeoAxis(
        f[1, 1];
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

    plottable = _convert_plottable(gdf, geom_col)
    plottable_buffer = _convert_plottable(gdf, buffer_col)

    gdf[:, color_by] .= ifelse.(ismissing.(gdf[:, color_by]), "NA", gdf[:, color_by])
    gdf[:, color_by] = convert.(String, gdf[:, color_by])

    # Define the unique colors and names for each level of factor color_by.
    # Use a different color palette for factors with high numbers of levels
    # (this palette is not as good for visualisation).
    if size(unique(gdf[:, color_by]),1) <= 20
        palette = ColorSchemes.tableau_20.colors
    else
        palette = ColorSchemes.flag_ec.colors
    end

    color_indices = groupindices(groupby(gdf, color_by))
    names = unique(DataFrame(indices=color_indices, names=gdf[:, color_by]))

    # Create the unique legend entries for each level of color_by
    unique_names = names.names
    legend_entries = []
    for name in eachrow(names)
        col = palette[name.indices]
        LE = PolyElement(; color=col)
        push!(legend_entries, [LE])
    end

    poly!(ga, plottable_buffer, color=palette[color_indices])
    scatter!(ga, plottable, transparency=true)

    Legend(f[2, 1], legend_entries, unique_names, nbanks=3, tellheight=true,
    tellwidth=false, orientation=:horizontal, labelsize=10)

    display(f)
end

plot_map(relevant_points; geom_col=:geometry, color_by=:facility, buffer_col=:buffer)

# Save the GBR-relevant points and the buffer zones for points. This crs is not the same as
# MPA crs in all regions. May need reprojection before visualisation/use with each region
GDF.write(joinpath(MPA_DATA_DIR, "../updated_boating_facilities.gpkg"), relevant_points, crs=GFT.EPSG(4326))
