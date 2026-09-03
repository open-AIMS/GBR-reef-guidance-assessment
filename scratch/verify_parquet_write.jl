"""
Standalone execution check for the new Parquet write + round-trip verification code
added to `valid_lookup` (`src/geom_handlers/lookup_processing.jl`), Phase 1 of
`.claude/plans/duckdb-switch.md`.

Loads the real, already-generated `Cairns-Cooktown_valid_slopes_lookup.arrow` file as
`area_store` (avoiding an expensive from-scratch raster recompute, which is outside this
track's scope) and runs the exact same write+verify logic `valid_lookup` runs, against
the real production data.
"""

using DataFrames
import Arrow
import QuackIO

const DATA_DIR = "/home/takuyai/projects/ReefGuide/ReefGuideWorker.jl/data"
const dst_file = joinpath(DATA_DIR, "Cairns-Cooktown_valid_slopes_lookup.arrow")
const parquet_dst_file = replace(
    dst_file, "_valid_slopes_lookup.arrow" => "_valid_slopes_lookup.parquet"
)

println("Loading real area_store from $(dst_file) ...")
@time area_store = DataFrame(Arrow.Table(dst_file))
println("Loaded: $(nrow(area_store)) rows, columns: $(names(area_store))")
println("Column types:")
for col in names(area_store)
    println("  $(col) => $(eltype(area_store[!, col]))")
end

println("\nWriting Parquet copy to $(parquet_dst_file) ...")
@time QuackIO.write_table(parquet_dst_file, area_store; format=:parquet)

println("\nRunning round-trip verification ...")
@time let
    roundtrip = QuackIO.read_parquet(DataFrame, parquet_dst_file)

    if nrow(roundtrip) != nrow(area_store)
        error(
            "Parquet round-trip verification failed for $(parquet_dst_file): " *
            "row count $(nrow(roundtrip)) != $(nrow(area_store))"
        )
    end
    println("  row count OK: $(nrow(roundtrip)) == $(nrow(area_store))")

    orig_cols = Set(names(area_store))
    rt_cols = Set(names(roundtrip))
    if orig_cols != rt_cols
        error(
            "Parquet round-trip verification failed for $(parquet_dst_file): " *
            "column set mismatch (original: $(orig_cols), round-trip: $(rt_cols))"
        )
    end
    println("  column set OK: $(orig_cols)")

    bounds_skip_cols = Set([:lon_idx, :lat_idx, :lons, :lats])
    for col in names(area_store)
        Symbol(col) in bounds_skip_cols && continue
        orig_vals = area_store[!, col]
        nonmissingtype(eltype(orig_vals)) <: Real || continue
        orig_non_missing = collect(skipmissing(orig_vals))
        isempty(orig_non_missing) && continue

        rt_non_missing = collect(skipmissing(roundtrip[!, col]))
        if minimum(orig_non_missing) != minimum(rt_non_missing) ||
            maximum(orig_non_missing) != maximum(rt_non_missing)
            error(
                "Parquet round-trip verification failed for $(parquet_dst_file): " *
                "min/max mismatch for column $(col)"
            )
        end
        println(
            "  $(col) min/max OK: ($(minimum(orig_non_missing)), $(maximum(orig_non_missing)))"
        )
    end

    n_spot_check = min(100, nrow(area_store))
    spot_rows = rand(1:nrow(area_store), n_spot_check)
    for row in spot_rows
        for col in names(area_store)
            if !isequal(area_store[row, col], roundtrip[row, col])
                error(
                    "Parquet round-trip verification failed for $(parquet_dst_file): " *
                    "value mismatch at row $(row), column $(col) " *
                    "($(area_store[row, col]) != $(roundtrip[row, col]))"
                )
            end
        end
    end
    println("  spot-check OK: $(n_spot_check) random rows compared equal (isequal, incl. missing)")

    roundtrip = nothing
end

println("\nAll verification checks passed.")
