using Statistics

using FLoops

import ArchGDAL as AG
import GeoInterface as GI
import GeometryOps as GO
using LibGEOS

using CoordinateTransformations

using Rasters


include("geom_ops.jl")

function identify_potential_sites(raster, threshold, geom, buffer_area)
    res = abs.(step.(dims(raster, (X, Y))))

    geom_buff = GO.buffer(geom, 10)

    potential_loc_idx = findall(raster .> threshold)
    raster = rebuild(raster; missingval=0)

    # Define maximum extent to consider
    lat1, lon1 = Int64.(Tuple(minimum(potential_loc_idx)) .- 1)
    lat2, lon2 = Int64.(Tuple(maximum(potential_loc_idx)) .+ 1)
    lon1 = dims(raster, X)[lon1]
    lon2 = dims(raster, X)[lon2]
    lat1 = dims(raster, Y)[lat1]
    lat2 = dims(raster, Y)[lat2]

    if lon1 > lon2
        swap = lon1
        lon1 = lon2
        lon2 = swap
    end

    if lat1 > lat2
        swap = lat1
        lat1 = lat2
        lat2 = swap
    end

    # Crop raster down to assessment area
    # TODO: Is this entire step just `trim(mask())` with missingval set to false/zero?
    assess_raster = view(raster, X(lon1 .. lon2), Y(lat1 .. lat2))
    @assert all(size(assess_raster) .> (0, 0))

    best_score = zeros(length(potential_loc_idx))
    best_poly = Vector(undef, length(potential_loc_idx))
    best_degree = zeros(Int64, length(potential_loc_idx))
    for (i, loc_idx) in enumerate(potential_loc_idx)
        # Create rectangle
        lon_lat = Tuple(loc_idx)
        lat = dims(raster, Y)[lon_lat[1]]
        lon = dims(raster, X)[lon_lat[2]]

        # Move geom to target centroid
        tf_geom = move_geom(geom_buff, (lon, lat))
        b_score, b_deg, b_poly = assess_reef_site(assess_raster, tf_geom)

        best_score[i] = b_score
        best_degree[i] = b_deg
        best_poly[i] = b_poly
    end

    return best_score, best_degree, best_poly
end


function assess_reef_site(rst, geom; degree_step=15.0)
    n_rotations = length(0.0:degree_step:359.0)
    score = zeros(n_rotations)
    best_poly = Vector(undef, n_rotations)

    # f, ax, sp = heatmap(rst)
    for (j, r) in enumerate(0.0:degree_step:359.0)
        rot_geom = rotate_geom(geom, r)
        c_rst = crop(rst; to=rot_geom)
        if !all(size(c_rst) .> (0, 0))
            @warn "No data found!"
            continue
        end

        score[j] = mean(c_rst)
        best_poly[j] = rot_geom
        # poly!(rot_geom, alpha=0.2)
    end

    # ax.title = "Polygon rotation"
    # resize_to_layout!(f)
    # display(f)

    return score[argmax(score)], argmax(score)-1, best_poly[argmax(score)]
end