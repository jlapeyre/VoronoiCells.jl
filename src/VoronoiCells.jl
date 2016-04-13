module VoronoiCells

# We do not compile
# 1. VoronoiDelaunay is not compiled.
# 2. VoronoiCells itself loads in about 100ms

include("voronoi_cells.jl")
include("store_cells.jl")

end # module
