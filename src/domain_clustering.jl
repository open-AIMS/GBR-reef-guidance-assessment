import GeoDataFrames as GDF
import ArchGDAL as AG
using DataFrames, StatsBase, Statistics
using Clustering, Distances
using BlackBoxOptim, Random


"""Identify column IDs which hold Geometry data types."""
function get_geometry_col(df::DataFrame)
    col_types = string.(eltype.(eachcol(df)))
    return first(findall(contains.(col_types, "Geometry")))
end

"""Get the geometry data as a vector"""
function get_geometry(df::DataFrame)::Vector
    geom_col = get_geometry_col(df)
    if sum(geom_col) == 0
        error("No geometry data found")
    end

    return df[:, geom_col]
end


"""
    centroids(df::DataFrame)

Extract and return long/lat from a GeoDataFrame.

# Arguments
- `df` : GeoDataFrame

# Returns
Array of tuples (x, y), where x and y relate to long and lat respectively.
"""
function centroids(df::DataFrame)::Vector{Tuple{Float64,Float64}}
    site_centroids::Vector = AG.centroid.(get_geometry(df))
    return collect(zip(AG.getx.(site_centroids, 0), AG.gety.(site_centroids, 0)))
end


"""
    create_distance_matrix(site_data::DataFrame)::Matrix

Calculate matrix of unique distances between locations.

# Returns
Distance between locations in meters
"""
function create_distance_matrix(site_data::DataFrame)::Matrix{Float64}
    site_centroids = centroids(site_data)
    longitudes = first.(site_centroids)
    latitudes = last.(site_centroids)

    n_sites = size(site_data, 1)
    dist = zeros(n_sites, n_sites)
    for ii in axes(dist, 1)
        for jj in axes(dist, 2)
            if ii == jj
                continue
            end

            @views dist[ii, jj] = haversine(
                (longitudes[ii], latitudes[ii]), (longitudes[jj], latitudes[jj])
            )
        end
    end

    return dist
end

"""Generate DxN point data based on location centroids."""
function create_location_points(gdf::DataFrame)
    all_centroids = centroids(gdf)
    all_locs = Matrix([first.(all_centroids) last.(all_centroids)]')

    return all_locs
end


function opt_kmeans(X, locations, dists)
    n_clusters = X[1]

    local clusters
    try
        clusters = kmeans(locations, floor(Int64, n_clusters), display=:none)
    catch err
        if err isa BoundsError
            return 1.0
        end
    end

    if !clusters.converged
        # Return worst score if k-means has not converged
        return 1.0
    end

    assignments = clusters.assignments

    sil_score = -1.0
    try
        sil_score = silhouettes(assignments, dists)
    catch err
        if !(err isa ArgumentError)
            rethrow(err)
        end

        # All locations assigned to a single cluster so assign worst score
        sil_score = -1.0
    end

    # Optimization direction is toward the minimum, so invert score.
    return -median(sil_score)
end

"""
    cluster(gpkg_path::String, output_path::String)

Cluster locations indicated in the geopackage file.

The floor of the

# Arguments
- `gpkg_path` : path to geopackage
- `output_path` : desired path to output geopackage
"""
function cluster(gpkg_path::String, output_path::String)
    Random.seed!(101)

    gdf = GDF.read(gpkg_path)
    dist_mat = create_distance_matrix(gdf)
    all_locs = create_location_points(gdf)

    n_locs = length(unique(gdf.UNIQUE_ID))
    dist_bnds = [
        # n clusters to find (n_locs + 1 to 3*n_reefs)
        # Upper bound has +1 as we take the floor of the sampled value.
        (n_locs+1, (n_locs*3)+1),
    ]

    opt_func = x -> opt_kmeans(x, all_locs, dist_mat)

    res = bboptimize(
        opt_func;
        SearchRange=dist_bnds,
        MaxSteps=100_000,
        ϵ=0.2,
    )
    best_params = best_candidate(res)

    clusters = kmeans(all_locs, floor(Int64, best_params[1]), display=:none)
    assignments = clusters.assignments
    sil_score = silhouettes(assignments, dist_mat)
    @info "Produced clusters with median Silhouette score of: $(median(sil_score))"
    @info "Silhouette score ranges from -1 to 1, with higher values being more desirable."
    @info "Values above 0.7 are regarded as strong, values over 0.5 are reasonable, and"
    @info "values below 0.25 are weak."
    #Qualitative indicators taken from https://en.wikipedia.org/wiki/Silhouette_(clustering)

    gdf.cluster_id = assignments
    GDF.write(output_path, gdf; geom_columns=(:geom, ))
end

cluster(joinpath(MPA_QGIS_DIR, "reef_suitability.gpkg"), joinpath(MPA_OUTPUT_DIR, "clustered_reef_suitability.gpkg"))
