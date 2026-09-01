# Standalone verification for the nullable-column behaviour added to `stack_values`
# (commit "Bathymetry-bound valid-pixel mask, nullable per-layer values").
#
# `stack_values` now allocates `Vector{Union{Missing,T}}` per layer and writes
# `missing` wherever a sampled value `isequal` that layer's own `missingval`
# sentinel, so a pixel valid under the bathymetry-bound mask can still be recorded
# as `missing` for a layer that does not cover it. This builds a tiny synthetic
# `RasterStack` with distinct per-layer sentinels (numeric, class, NaN, all-missing)
# and checks that behaviour, the row-scatter ordering, and the Arrow round-trip,
# without needing region data or `common.jl`'s `.config.toml`.
#
# Not covered: the band-splitting loop runs once here (the `band_height % 256 == 0`
# assertion rules out a multi-band fixture this small).
#
# Run:  julia --project=. scratch/verify_missing_values.jl

using Test
using SparseArrays
using Rasters
using Arrow: Arrow
using DataFrames

# `lookup_processing.jl` calls `force_gc_cleanup()` (defined in `common.jl`); stub it
# so the file can be included on its own.
force_gc_cleanup()::Nothing = nothing

include(joinpath(@__DIR__, "..", "src", "geom_handlers", "lookup_processing.jl"))

# 5 (X) by 4 (Y) grid on a regular 10-unit step; centroids are cell-corner + 5.
xs = 0.0:10.0:40.0
ys = 0.0:10.0:30.0
xdim = X(xs)
ydim = Y(ys)

const NUM_SENTINEL = -9999.0f0
const CLASS_SENTINEL = Int8(0)

# Depth is fully covered; depth[x, y] == (y - 1) * 5 + x by column-major reshape.
depth = Raster(Float32.(reshape(1:20, 5, 4)), (xdim, ydim); missingval=NUM_SENTINEL)

turbidity = Raster(fill(2.5f0, 5, 4), (xdim, ydim); missingval=NUM_SENTINEL)
turbidity[2, 3] = NUM_SENTINEL          # a valid pixel this layer does not cover

benthic = Raster(fill(Int8(7), 5, 4), (xdim, ydim); missingval=CLASS_SENTINEL)
benthic[1, 2] = CLASS_SENTINEL          # class layer uses its own sentinel

waves = Raster(fill(1.0f0, 5, 4), (xdim, ydim); missingval=NaN32)
waves[3, 1] = NaN32                     # NaN sentinel: relies on `isequal`, not `==`

all_missing = Raster(fill(NUM_SENTINEL, 5, 4), (xdim, ydim); missingval=NUM_SENTINEL)

rst_stack = RasterStack((;
    Depth=depth, Turbidity=turbidity, Benthic=benthic, Waves=waves, AllMissing=all_missing
))

# Valid pixels as (x_idx, y_idx). The mask is walked Y-major (CSC order) but rows are
# emitted sorted by (x_idx, y_idx), so the scatter permutation here is
# [1, 4, 2, 3, 5] - a 3-cycle, not an involution, so a `sortperm` used where
# `invperm(sortperm(...))` is meant would reorder the output and fail below.
# (3, 1) and (3, 3) share an X, exercising the Y tiebreaker in the sort.
valid_pixels = [(1, 1), (3, 1), (1, 2), (2, 3), (3, 3)]
valid_mask = sparse(
    first.(valid_pixels), last.(valid_pixels), trues(length(valid_pixels)), 5, 4
)

# Emitted row order: (x_idx, y_idx) ascending.
emitted = sort(valid_pixels)

value_cols, lon_idx_col, lat_idx_col, centroid_lons, centroid_lats = stack_values(
    valid_mask, rst_stack
)

@testset "stack_values nullable per-layer columns" begin
    @testset "column eltypes are nullable" begin
        @test eltype(value_cols.Depth) == Union{Missing,Float32}
        @test eltype(value_cols.Turbidity) == Union{Missing,Float32}
        @test eltype(value_cols.Benthic) == Union{Missing,Int8}
        @test eltype(value_cols.Waves) == Union{Missing,Float32}
        @test eltype(value_cols.AllMissing) == Union{Missing,Float32}
    end

    @testset "rows are emitted sorted by (x_idx, y_idx)" begin
        @test collect(zip(lon_idx_col, lat_idx_col)) == emitted
    end

    @testset "centroids are cell-corner + half a pixel, in emitted row order" begin
        @test centroid_lons == [xs[x] + 5.0 for (x, _) in emitted]
        @test centroid_lats == [ys[y] + 5.0 for (_, y) in emitted]
    end

    @testset "sentinel values become missing, everything else is verbatim" begin
        # Expected columns in emitted row order.
        @test isequal(value_cols.Depth, Float32[(y - 1) * 5 + x for (x, y) in emitted])
        @test isequal(
            value_cols.Turbidity, [(x, y) == (2, 3) ? missing : 2.5f0 for (x, y) in emitted]
        )
        @test isequal(
            value_cols.Benthic, [(x, y) == (1, 2) ? missing : Int8(7) for (x, y) in emitted]
        )
        @test isequal(
            value_cols.Waves, [(x, y) == (3, 1) ? missing : 1.0f0 for (x, y) in emitted]
        )
        @test all(ismissing, value_cols.AllMissing)
    end

    @testset "Arrow round-trip preserves nullability" begin
        df = DataFrame(; lon_idx=lon_idx_col, lat_idx=lat_idx_col, value_cols...)
        path = tempname() * ".arrow"
        Arrow.write(path, df; compress=:zstd)
        back = DataFrame(Arrow.Table(path))
        rm(path)

        @test eltype(back.Turbidity) == Union{Missing,Float32}
        @test eltype(back.Benthic) == Union{Missing,Int8}
        @test isequal(back.Turbidity, value_cols.Turbidity)
        @test isequal(back.Benthic, value_cols.Benthic)
        @test isequal(back.Waves, value_cols.Waves)
        @test all(ismissing, back.AllMissing)
        @test count(ismissing, back.Depth) == 0
    end
end
