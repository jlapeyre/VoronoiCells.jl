# VoronoiCells

This module provides access to Vornoi cells in the tesselation created
by [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl)

This the first version, with a rather sprawling interface that should not be
considered stable. Following is a brief description of some of the
most recent, high-level interface.

### Some functions

```julia
cells = poissonvoronoicells(npts)
```
Sample (approximately) a Poisson point process with exactly `npts` points.
Return the Voronoi cells. Only complete cells within a square area are returned.
The cells very near the boundary do not have the same geometrical statistics
as those in the bulk.

Lookup of the cell containing a point is done via a square grid. By default, about
100 generator points are put in each square. `cells = poissonvoronoicells(npts,ndiv)`,
will put about `ndiv^2` generator points in each square, at the cost of a higher storage
requirement with increasing `ndiv`.

```julia
smaxcoord(cells), smincoord(cells)
```
Scaled maximum and minimum of both `x` and `y` coordinates. The scaling is such that
the average cell size is (approximately) `1` and the origin is at the center of the square area.

```julia
idx = sfindindex(cells,x,y)
```
Search efficiently for and return the index of the cell containing the scaled point `(x,y)`.
The cell may then be retrieved via `cells[idx]`.

```julia
isvalid(idx)
```
Return `true` if a cell was found containing `(x,y)`, otherwise `false`.

```julia
c = cells[idx], sarea(cells,idx), nedges(c), nverts(c)
```
Scaled area of cell at index `idx`. Number of edges and number of vertices in the cell.

```julia
cells = poissonvoronoicells(npts)
write("fname.dat", cells)
```
Write cells to a file. The file is in a fast binary format. It does not
use any generic serialization.

```julia
cells = read_one_cellfile("fname.dat")
```
Read cells from file written by `write`.

### Example

Here is sample code that inefficiently measures the second moment of the cell area
distribution by doing a random walk and recording the area of the containing cell
at each step.
```julia
flip() = 2*rand(Bool) - 1

function walkmean()
    npts = 100000  # number of points in Poisson point process
    cells = poissonvoronoicells(npts)  # cells with unit mean size
    nsteps = 10^6  # number of walk steps
    dx = 0.1       # step length
    x = 0.0        # walker's position
    y = 0.0
    sumarea = 0.0  # sum of areas of cells at each step
    nareas = 0     # number of area samples collected
    for stp in 1:nsteps
        x += flip() * dx
        y += flip() * dx
        idx = sfindindex(cells,x,y)  # find index of cell containing point (x,y)
        if ! isvalid(idx)
            println("Exited tesselated area at step number $stp.")
            break
        end
        sumarea += sarea(cells,idx)  # add (scaled) area of cell
        nareas += 1
    end
    @printf("mean sampled area %.4f,  nareas %d\n", sumarea/nareas,  nareas)
end
```

There is more to the interface including unscaled versions of functions. The
unscaled area is a square with both `x` and `y` coordinates restricted
to approximately `1 <= x,y <= 2`. See [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl)
for details. More of the interface may be documented later.

```julia
sizeof(poissonvoronoicells(npts))
```

Size in bytes for storage of cells. This includes cell generators, edges, and the grid structure.

```julia
cells = poissonvoronoicells(npts)
sizeof(cells[1])
```
sizeof a single cell.
