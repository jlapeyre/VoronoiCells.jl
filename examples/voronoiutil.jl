using VoronoiDelaunay
using VoronoiCells

# does not work due to bitrot
function printpoissonvoronoicells(fname,n)
    cells = generate_voronoi_cells(n)
    print_cells(fname,cells)
    cells
end
