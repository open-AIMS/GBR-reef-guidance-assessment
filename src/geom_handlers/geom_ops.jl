using Statistics

using FLoops

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO

using CoordinateTransformations, Rotations

using Rasters

using GLMakie


function create_poly(verts, crs)
    sel_lines = GI.LineString(GI.Point.(verts))
    ring = GI.LinearRing(GI.getpoint(sel_lines))

    return GI.Polygon([ring]; crs=crs)
end

"""
    create_bbox(xs::Tuple, ys::Tuple)::Vector{Tuple{Float64, Float64}}

Create bounding box from x and y coordinates

Returns in order of top left, top right, bottom right, bottom left
"""
function create_bbox(xs::Tuple, ys::Tuple)::Vector{Tuple{Float64,Float64}}
    # Top left, top right, bottom right, bottom left
    return [
        (xs[1], ys[2]),
        (xs[2], ys[2]),
        (xs[2], ys[1]),
        (xs[1], ys[1]),
        (xs[1], ys[2])
    ]
end

"""Rotate the polygon by the given angle about its center."""
function rotate_polygon(poly_points, centroid, degrees)
    if degrees == 0.0
        return poly_points
    end

    theta = deg2rad(degrees)
    cosang, sinang = cos(theta), sin(theta)

    # Center is used as pivot point
    cx, cy = centroid

    # Update the coordinates of each vertex
    new_points = copy(poly_points)
    for (i, p) in enumerate(poly_points)
        x, y = p
        x -= cx
        y -= cy
        new_x = x * cosang - y * sinang + cx
        new_y = x * sinang + y * cosang + cy

        new_points[i] = (new_x, new_y)
    end

    return new_points
end

function rotate_geom(geom, degrees::Float64)
    if degrees == 0.0
        return geom
    end

    theta = deg2rad(degrees)
    cosang, sinang = cos(theta), sin(theta)

    # Center is used as pivot point
    cx, cy = GO.centroid(geom)

    # Extract points
    new_points = try
        getfield.(GI.getpoint(geom), :geom)
    catch err
        if !contains(err.msg, "type SArray has no field geom")
            throw(err)
        end

        collect(GI.getpoint(geom))
    end

    # Calculate new coordinates of each vertex
    for (i, p) in enumerate(new_points)
        x, y = p
        x -= cx
        y -= cy
        new_x = x * cosang - y * sinang + cx
        new_y = x * sinang + y * cosang + cy

        new_points[i] = (new_x, new_y)
    end

    return create_poly(new_points, GI.crs(geom))
end

function move_geom(geom, new_centroid::Tuple)
    tf_lon, tf_lat = new_centroid .- GO.centroid(geom)
    f = CoordinateTransformations.Translation(tf_lon, tf_lat)
    return GO.transform(f, geom)
end


# using GLMakie, GeoMakie

# include("geom_handlers/site_assessment.jl")


# rst = Raster("../outputs/MPA/_Cairns-Cooktown_suitable_flats_new.tif")
# Flip so that longs are along the X dimension, and lats are along the Y dimension
# Note: This is largely unnecessary, and was really on to make X = lon and Y = lat
# rst2 = resample(rst; crs=EPSG(7856))'

# # b_score, b_degree, b_polys = identify_potential_sites(rst2, 80.0, 10, 150, EPSG(7856))

# # Define the polygon shape to search for (and auto-rotate)
# xs = (1, 450)
# ys = (1, 10)
# search_plot = create_poly(create_bbox(xs, ys), EPSG(7856))

# b_score, b_degree, b_polys = identify_potential_sites(rst2, 80.0, search_plot, 5.0)
