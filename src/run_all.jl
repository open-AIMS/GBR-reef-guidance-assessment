ENV["JULIA_DEBUG"] = "Main"  # Enable debug messages

include("01a_prep_hybrid_benthic.jl")
include("01b_prep_hybrid_geomorphic.jl")
include("02_prep_MPA.jl")

@time begin
    @time prep_hybrid_benthic()
    @time prep_hybrid_geomorphic()
    @time prep_MPA()
end
