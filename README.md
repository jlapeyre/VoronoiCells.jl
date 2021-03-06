# VoronoiCells
*provides access to Voronoi cells in [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl)*

Linux, OSX: [![Build Status](https://travis-ci.org/jlapeyre/VoronoiCells.jl.svg)](https://travis-ci.org/jlapeyre/VoronoiCells.jl)
&nbsp;
Windows: [![Build Status](https://ci.appveyor.com/api/projects/status/github/jlapeyre/VoronoiCells.jl?branch=master&svg=true)](https://ci.appveyor.com/project/jlapeyre/voronoicells-jl)
&nbsp; &nbsp; &nbsp;
[![Coverage Status](https://coveralls.io/repos/jlapeyre/VoronoiCells.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/jlapeyre/VoronoiCells.jl?branch=master)
[![codecov.io](http://codecov.io/github/jlapeyre/VoronoiCells.jl/coverage.svg?branch=master)](http://codecov.io/github/jlapeyre/VoronoiCells.jl?branch=master)

#### Other package with same name

There is a package by the same name [VoronoiCells.jl](https://github.com/JuliaGeometry/VoronoiCells.jl) that was created about
nine months later than this package. The other package is registered, this one is not.

#### Installation

`VoronoiCells` depends on `VoronoiDelaunay.jl` which depends on `GeometricalPredicates.jl`, and often, only the
development version of `GeometricalPredicates.jl` will successfully load. You can switch to the development version like this:

```julia
Pkg.add("GeometricalPredicates")
Pkg.checkout("GeometricalPredicates")
Pkg.clone("https://github.com/jlapeyre/VoronoiCells.jl")
```

### Some functions

```julia
cells = poissonvoronoicells(npts)
```
Sample a Poisson point process in the plane with mean `npts` 'generator' points per unit square and
return the Voronoi cells. Only complete cells within the unit square [1,2] x [1,2] are returned.
The cells very near the boundary do not have the same geometrical statistics
as those in the bulk. `length(cells)` returns the number of cells retained. If `Distributions.jl` is
installed, then the number of generators in the square is sampled from the Poisson distribution, otherwise,
the number is taken to be exactly `npts`. Use  `approxpoissonvoronoicells(npts)` to force the latter behavior.

Lookup of the cell containing a point is done via a square grid. By default, about
100 generator points are put in each square. `cells = poissonvoronoicells(npts,ndiv)`,
will put about `ndiv^2` generator points in each square, at the cost of a higher storage
requirement with increasing `ndiv`.

```julia
cells = poissonvoronoicellsnogrid(npts)
```

Return an array of `VoronoiCell`'s without the grid structure for locating them.


```julia
smaxcoord(cells), smincoord(cells)
```
Scaled maximum and minimum of both `x` and `y` coordinates. The scaling is such that
the average cell size is (approximately) `1` and the origin is at the center of the square area.


search efficiently for and returns the index of the cell containing the scaled point `(x,y)`.
The cell may then be retrieved via `cells[idx]`.  `findindex` finds cells at unscaled (within the shifted unit square)
coordinates. If no cell contains the point `(x,y)`, an invalid index is returned. This condition is checked with `isvalid()`.

```julia
isvalid(idx)
```
Return `true` if a cell was found containing `(x,y)`, otherwise `false`.

```julia
getinvalidcellindex()
```
Return an invalid cell index.

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

```julia
cells = poissonvoronoicells(npts)
isexternal(cell[i])
```
Return `true` if the vertices of the cell lie entirely within the unit square [1,2] x [1,2].


```julia
cells = poissonvoronoicells(npts)
random_cell(cell)
```
Return the cell at a randomly chosen point on the plane. (Resampling until a point
containing a cell is found.

```julia
t = poissontesselation(npts)
t = approxpoissontesselation(npts)
```
Return Delaunay tesselation of the Poisson point process in the plane.

```julia
gcells = voronoicells(t)
cells = voronoicellsnogrid(t)
```
Return the cells corresponding to the `DelaunayTessellation2D` `t`.

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

### Unexported functions

These functions have doc strings.


<!--  LocalWords:  VoronoiCells Voronoi VoronoiDelaunay jl OSX nbsp
 -->
<!--  LocalWords:  codecov io GeometricalPredicates julia npts ndiv
 -->
<!--  LocalWords:  poissonvoronoicells approxpoissonvoronoicells idx
 -->
<!--  LocalWords:  poissonvoronoicellsnogrid VoronoiCell smaxcoord dx
 -->
<!--  LocalWords:  smincoord sfindindex findindex isvalid searche dat
 -->
<!--  LocalWords:  unscaled getinvalidcellindex sarea nedges nverts
 -->
<!--  LocalWords:  fname cellfile isexternal Resampling Delaunay Bool
 -->
<!--  LocalWords:  poissontesselation approxpoissontesselation gcells
 -->
<!--  LocalWords:  tesselation voronoicells voronoicellsnogrid nsteps
 -->
<!--  LocalWords:  DelaunayTessellation walkmean sumarea nareas stp
 -->
<!--  LocalWords:  println tesselated printf sizeof
 -->
