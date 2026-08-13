"""
Audit deployed `*_valid_slopes_lookup.arrow` files for the pre-computed `:lons`/`:lats`
columns added by `valid_lookup()` (see `geom_handlers/lookup_processing.jl`).

Any lookup file written before that column was added is stale: `valid_lookup()` skips
regeneration whenever its `dst_file` already exists, so a stale file will silently keep
forcing `ReefGuide.jl` down the slower runtime lat/lon fallback path
(`add_lat_long_columns_to_dataframe`) forever. This script finds and regenerates them.

Usage: `julia --project src/audit_slope_lookups.jl` from the package root, after `.config.toml`
is set up per the README.
"""

include("02_prep_MPA.jl")  # pulls in common.jl and criteria path/lookup helpers

"""
    audit_and_regenerate_slope_lookups(; regenerate::Bool=true)::Nothing

Scan `MPA_OUTPUT_DIR` for `*_valid_slopes_lookup.arrow` files, report which are missing the
`:lons`/`:lats` columns, and (when `regenerate=true`) delete and rebuild each stale file via
`valid_lookup()`.

# Arguments
- `regenerate` : If `false`, only report stale files without deleting/rebuilding them.
"""
function audit_and_regenerate_slope_lookups(; regenerate::Bool=true)::Nothing
    lookup_files = glob("*_valid_slopes_lookup.arrow", MPA_OUTPUT_DIR)
    if isempty(lookup_files)
        @info "No *_valid_slopes_lookup.arrow files found under $(MPA_OUTPUT_DIR)."
        return nothing
    end

    stale_files = String[]
    for fn in lookup_files
        area_store = DataFrame(Arrow.Table(fn))
        if !("lons" in names(area_store) && "lats" in names(area_store))
            push!(stale_files, fn)
        end
    end

    if isempty(stale_files)
        @info "All $(length(lookup_files)) slope lookup file(s) already have :lons/:lats columns."
        return nothing
    end

    @warn "$(length(stale_files)) of $(length(lookup_files)) slope lookup file(s) are missing :lons/:lats and will fall back to the slower runtime lat/lon computation:" stale_files

    if !regenerate
        return nothing
    end

    for fn in stale_files
        # Recover the region name from the standard "<region>_valid_slopes_lookup.arrow" naming.
        reg = replace(basename(fn), "_valid_slopes_lookup.arrow" => "")
        if reg ∉ REGIONS
            @warn "Could not map $(fn) back to a known region (parsed as $(reg)); skipping."
            continue
        end

        @info "Regenerating $(fn) for region $(reg)..."
        criteria_paths = create_criteria_paths(reg)
        valid_slopes_fn = joinpath(MPA_OUTPUT_DIR, "$(reg)_valid_slopes.tif")

        rm(fn)
        valid_lookup(criteria_paths, valid_slopes_fn, fn)
    end

    return nothing
end

audit_and_regenerate_slope_lookups()
