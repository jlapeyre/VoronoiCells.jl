using VoronoiDelaunay
using VoronoiCells

function poissontesselation(n)
    width = max_coord - min_coord
    a = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:n]
    tess = DelaunayTessellation()
    push!(tess,a)
    tess
end

function poissonvoronoicells(n,ngrid)
    tess = generate_poisson_tesselation(n)
    voronoicells(tess,ngrid)
end

function poissonvoronoicells(n)
    ngrid = round(sqrt(n)/10)
    poissonvoronoicells(n,ngrid)
end

function printpoissonvoronoicells(fname,n)
    cells = generate_voronoi_cells(n)
    print_cells(fname,cells)
    cells
end

