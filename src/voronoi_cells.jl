using VoronoiDelaunay
using GeometricalPredicates
import GeometricalPredicates
import GeometricalPredicates: area  # needed so that area() can be called in runtests.jl
import VoronoiDelaunay.isexternal, VoronoiDelaunay.locate, VoronoiDelaunay.findindex
import Base: getindex, ==, -, sizeof, endof

if Pkg.installed("Distributions") == nothing
    @eval const haveDistributions = false
else
    @eval const haveDistributions = true
    using Distributions
end


export VoronoiCell, VoronoiCellsA, voronoicells, voronoicellsnogrid, findindex, locate, invoronoicell, area, avoronoicellsnogrid
export getcellindex,isexternal, nverts, nedges, getgenerator, scale, iscale, get_

export VoronoiCellIndex, isvalid, getinvalidcellindex

# For VoronoiCellsA
export cellstogrid, poissontesselation, poissonvoronoicells, poissonvoronoicellsnogrid,
approxpoissontesselation, approxpoissonvoronoicells, approxpoissonvoronoicellsnogrid,
       ngrid, findindex0
export slocate, sfindindex0, sfindindex, area, sarea, smaxcoord, smincoord, sgetgenerator
export getareascale, getscale, getshift, getcells, random_cell, random_cell_index

####

const min_coord = GeometricalPredicates.min_coord + eps(Float64)
const max_coord = GeometricalPredicates.max_coord - eps(Float64)

#### Some functions that only operate on GeometricalPredicates

# Euclidean distance
distance(p1::Point2D,p2::Point2D) = (getx(p2)-getx(p1))^2 + (gety(p2)-gety(p1))^2

# Maybe faster than Euclidean, and is sufficient for use below
distanceL1(p1::Point2D,p2::Point2D) = abs(getx(p2)-getx(p1)) + abs(gety(p2)-gety(p1))

# Removed use of this. Maybe remove this, as well
-(p1::Point2D,p2::Point2D) = Point2D(getx(p1)-getx(p2),gety(p1)-gety(p2))

# sign shows which side of line between P(x1,y1) P(x2,y2) the point

lineside(p1,p2,p) = (gety(p)-gety(p1))*(getx(p2)-getx(p1)) - (getx(p)-getx(p1))*(gety(p2)-gety(p1))

# true if p is inside poly, ie. if lineside gives the same sign for all edges in Voronoi cell
"""
    inconvexpolygon(poly::Array{Point2D,1}, p0::Point2D)

return `true` if `p0` is inside the polygon `poly`.
"""
function inconvexpolygon(poly::Array{Point2D,1}, p0::Point2D)
    is_inconvexpolygon = true
    n = length(poly)
    p1 = poly[n]
    p2 = poly[1]
    s = lineside(p1,p2,p0)
    @inbounds for i in 1:length(poly)-1
        p1 = poly[i]
        p2 = poly[i+1]
        s1::Float64 = lineside(p1,p2,p0)
        if sign(s1) != sign(s)
            is_inconvexpolygon = false
            break
        end
    end
    is_inconvexpolygon
end

function random_coordinate()
    width = max_coord - min_coord
    min_coord+rand()*width
end

function random_coordinate(cutoff)
    width = (max_coord - min_coord) - 2*cutoff
    min_coord + cutoff + rand()*width
end


random_point() = Point2D(random_coordinate(),random_coordinate())
random_point(cutoff) = Point2D(random_coordinate(cutoff),random_coordinate(cutoff))

#### Some functions that only operate on VoronoiDelaunay objects.

# Get vertex by integer index (1,2,3) are (a,b,c)
function get_vertex(trig::VoronoiDelaunay.DelaunayTriangle, iv::Int)
    if iv == 1
        geta(trig)
    elseif iv == 2
        getb(trig)
    elseif iv == 3
        getc(trig)
    else
        error("Vertex index must be 1,2, or 3. Got $iv.")
    end
end

# Get the neigboring triangle sharing side opposite vertex a,b,c via indices 1,2,3.
# Note regarding DelaunayTriangle:
# tr._neighbour_a is the neighboring triangle sharing the edge in tr that is opposite point a
# tr._neighbour_b, tr._neighbour_c contain the vertex a, while tr._neighbour_a does not.
function get_neighbor(trigs,trig,iv::Int)
    if iv == 1
        trigs[trig._neighbour_a]
    elseif iv == 2
        trigs[trig._neighbour_b]
    elseif iv == 3
        trigs[trig._neighbour_c]
    else
        error("Neighbor index must be 1,2, or 3. Got $iv.")
    end
end

# Same as above, but return index of triangle in array of Delaunay triangles,
# rather than the triangle itself
function get_neighbor_index(trigs,trig,iv::Int)
    if iv == 1
        trig._neighbour_a
    elseif iv == 2
        trig._neighbour_b
    elseif iv == 3
        trig._neighbour_c
    else
        error("Neighbor index must be 1,2, or 3. Got $iv.")
    end
end

# Return 1,2, or 3 if one of the three vertices a,b, or c
# in triangle tr is equal to point p0. error otherwise.
function match_vertex(tr, p0)
    eps0 = eps(min_coord) # same as for max_coord
    pa = geta(tr)
    if distanceL1(pa,p0) <= eps0
        return 1
    else
        pb = getb(tr)
        if distanceL1(pb,p0) <= eps0
            return 2
        else
            pc = getc(tr)
            if distanceL1(pc,p0) <= eps0
                return 3
            end
        end
    end
    error("match_vertex: Given triangle does not contain given vertex")
end

# return complement of {v1,v2} in {1,2,3}
# i.e. if input is 3,2, we return 1.
function match_other_vertex(v1,v2)
    v1 == v2 &&  error("Can't find other vertex: v1=$v1 v2=$v2")
    if v1 == 1
        if v2 == 2
            return 3
        else
            return 2
        end
    elseif v1 == 2
        if v2 == 1
            return 3
        else
            return 1
        end
    elseif v1 == 3    # wastes time, but safer
        if v2 == 1
            return 2
        else
            return 1
        end
    else
        error("Can't find other vertex: v1=$v1 v2=$v2")
    end
end

# Find the next DelaunayTriangle triangle when building a VoronoiCell
# trigs -- array of DelaunayTriangles
# tr -- current triangle
# iv_gen -- index (recall (1,2,3) <--> (a,b,c)) of generator point in the triangle tr
#  for the cell we are building. Each triangle has the generator as a vertex, but
#  the index may vary with the triangles.
# iv_opp -- index of the vertex opposite the side that is shared with the
#  next triangle.
# iv_other -- index of remaining vertex
# return:
# tr2  -- next triangle in the cell
# itr2  -- index into the array of all DelaunayTriangle's of tr2
# iv_gen2 -- index (1,2,3) in tr2 of vertex coincident with the generator.
# iv_opp2 -- index in tr2 of vertex opposite the edge in tr2 shared with
#  the next-next triangle.
# iv_other2 -- remaining vertex in tr2
@inline function next_cell_triangle(trigs,tr,iv_gen::Int,iv_opp::Int, iv_other::Int)
    eps0 = eps(min_coord)
    tr2 = get_neighbor(trigs,tr,iv_opp)
    itr2 = get_neighbor_index(trigs,tr,iv_opp)
    iv_gen2 = match_vertex(tr2,get_vertex(tr,iv_gen))
    iv_opp2 = match_vertex(tr2,get_vertex(tr,iv_other))
    iv_other2 = match_other_vertex(iv_gen2, iv_opp2)
    return (tr2, itr2, iv_gen2, iv_opp2, iv_other2)
end

#### Begin code specific to VoronoiCells

"""
    VoronoiCell

The VoronoiCell type. It contains a generator point and an array of the vertices of
the cell. The field `_generator` is one of the points in the point process. There is exactly one such point in each
cell.
"""
immutable VoronoiCell
    _generator::Point2D
    _verts::Array{Point2D,1}
    # neigbhors ?
end

VoronoiCell(generator) = VoronoiCell(generator, Array(Point2D, 0))
nverts(c::VoronoiCell) = length(c._verts)
nedges(c::VoronoiCell) = length(c._verts) + 1
getgenerator(c::VoronoiCell) = c._generator

# Area of irregular polygon. Don't make use of convex property.
"""
    area(c::VoronoiCell)

return the area of `c::VoronoiCell`.
"""
function area(c::VoronoiCell)
    vs::Array{Point2D,1} = c._verts
    A::Float64 = 0.0
    n::Int = length(vs)
   @inbounds @simd for i in 1:n-1
        A += getx(vs[i])*gety(vs[i+1]) -  getx(vs[i+1])*gety(vs[i])
    end
    A += getx(vs[n])*gety(vs[1]) -  getx(vs[1])*gety(vs[n])
    0.5 * abs(A)
end

# is any vertex on boundary of cell outside of allowed range
function isexternal(c::VoronoiCell)
    vs::Array{Point2D,1} = c._verts
    n::Int = length(vs)
    found_external::Bool = false
    @inbounds for i in 1:n
        xc::Float64 = getx(vs[i])
        if xc < min_coord || xc > max_coord
            found_external = true
            break
        end
        yc::Float64 = gety(vs[i])
        if yc < min_coord || yc > max_coord
            found_external = true
            break
        end
    end
    found_external
end

invoronoicell(c::VoronoiCell, p::Point2D) = inconvexpolygon(c._verts,p)

# Find the index in an array of cells of the (first) cell containing point p,
# or zero if no cell contains p.
# Use linear search. We do not use this function on type VoronoiCellsA below.
# It is too expensive to search all cells. Searching for one point in each of 10^6
# cells takes on the order of 1 day cpu time.
# Below, there is a method for this function for VoronoiCellsA.
function findindex(cells::Array{VoronoiCell,1}, p::Point2D)
    ifound::Int = 0
    @inbounds for i in 1:length(cells)
        if invoronoicell(cells[i],p)
            ifound = i
            break
        end
    end
    ifound
end

# could do this with iterator, I suppose, but that constructs all cells every time
# return the cell rather than index. Do not check for error (index == 0 )
locate(cells::Array{VoronoiCell,1}, p::Point2D) = cells[findindex(cells,p)]


function find_cell(trigs,tr, iv_gen, iv_opp, visited)
    cell::VoronoiCell = VoronoiCell(get_vertex(tr, iv_gen))
    iv_other = match_other_vertex(iv_gen,iv_opp)
    tr0 = tr
    is_visited::Bool = false
    while true
        (tr2, itr2, iv_gen, iv_opp, iv_other) = next_cell_triangle(trigs, tr, iv_gen, iv_opp, iv_other)
        if visited[itr2] && tr2 != tr0
            is_visited = true
            break
        end
        push!(cell._verts, circumcenter(tr))
        if tr2 == tr0 break end
        tr = tr2
    end
    (is_visited,cell)
end

# Using  produce, collect is much slower than just returning an array with avoronoicellsnogrid

"""
    avoronoicellsnogrid(t::DelaunayTessellation2D)

return an iterator of the cells `a::Array{VoronoiCell,1}` in the tesselation `t`.
This function uses `Task` and is much slower than `avoronoicellsnogrid(t::DelaunayTessellation2D)`, which
simply returns an array.
"""
function voronoicellsnogrid(t::DelaunayTessellation2D)
    visited::Array{Bool,1} = zeros(Bool, t._last_trig_index)
    visited[1] = true
    function voronoicelliterator()
        j::Int = 0
	@inbounds for ix::Int in 2:t._last_trig_index
	    visited[ix] && continue
	    const tr = t._trigs[ix]
	    visited[ix] = true
            isexternal(tr) && continue
            for iv_gen in 1:3

                iv_opp = mod1(iv_gen+1,3) # pick one of the other two vertices
                ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
                visited[ix2] && continue

                iv_opp = mod1(iv_gen+2,3) # pick the remaining vertex
                ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
                visited[ix2] && continue

                (is_visited, cell::VoronoiCell) = find_cell(t._trigs, tr, iv_gen, iv_opp, visited)
                is_visited && continue
                j += 1
#                j > 10000 && return
                isexternal(cell) && continue # just checking each triangle does not seem to work.
                produce(cell)
            end
        end
    end
    Task(voronoicelliterator)
end

"""
    avoronoicellsnogrid(t::DelaunayTessellation2D)

return an array of the cells `a::Array{VoronoiCell,1}` in the tesselation `t`.
This function is faster than `voronoicellsnogrid(t::DelaunayTessellation2D)`, which
uses lightweight threads.
"""
function avoronoicellsnogrid(t::DelaunayTessellation2D)
    visited = zeros(Bool, t._last_trig_index)
    visited[1] = true
    cells = Array(VoronoiCell,0)
    j = 0
    @inbounds for ix in 2:t._last_trig_index
	visited[ix] && continue
	const tr = t._trigs[ix]
	visited[ix] = true
        isexternal(tr) && continue
        @inbounds for iv_gen in 1:3

            iv_opp = mod1(iv_gen+1,3) # pick one of the other two vertices
            ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
            visited[ix2] && continue

            iv_opp = mod1(iv_gen+2,3) # pick the remaining vertex
            ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
            visited[ix2] && continue

            (is_visited, cell) = find_cell(t._trigs, tr, iv_gen, iv_opp, visited)
            is_visited && continue
            j += 1
            #                j > 10000 && return
            isexternal(cell) && continue # just checking each triangle does not seem to work.
            push!(cells,cell)
        end
    end
    cells
end

#### VoronoiCellsA

# Voronoi cells arranged in efficiently searchable array
# structure is _ngrid by _ngrid array of arrays of indices into the linear array.
# elements in the 2d array correspond to square regions covering area of point process.
#  of all VoronoiCell's.
# _cells -- linear array of all cells (currently in (almost) no order)
# _grid -- 2d array of arrays of indices into _cells.
# _scale, _shift, _areascale. For convenience of user. See scaled versions of functions below.
immutable VoronoiCellsA
    _ngrid::Int
    _cells::Array{VoronoiCell,1}
    _grid::Array{Array{Int,1},2}
    _scale::Float64
    _shift::Float64
    _areascale::Float64
end

# index into grid and array structure VoronoiCellsA
immutable VoronoiCellIndex
    _ix::Int
    _iy::Int
    _ind::Int
end

# Pick random points on the shifted unit square and return the cell there.
# This weights cells by their areas.
function random_cell(c::VoronoiCellsA)
    ntries::Int = 0
    local idx::VoronoiCellIndex
    while true
        pt = random_point()
        ntries += 1
        idx = findindex0(c,pt)
        if isvalid(idx) break end
        if ntries > 1000 error("VoronoiCells: can't find a random cell") end
    end
    c[idx]
end

function random_cell_index(c::VoronoiCellsA)
    ntries::Int = 0
    local idx::VoronoiCellIndex
    while true
        pt = random_point()
        ntries += 1
        idx = findindex0(c,pt)
        if isvalid(idx) break end
        if ntries > 1000 error("VoronoiCells: can't find a random cell") end
    end
    idx
end


function random_cell(c::VoronoiCellsA, cutoff)
    ntries::Int = 0
    local idx::VoronoiCellIndex
    while true
        pt = random_point(cutoff)
        ntries += 1
        idx = findindex0(c,pt)
        if isvalid(idx) break end
        if ntries > 1000 error("VoronoiCells: can't find a random cell") end
    end
    c[idx]
end

function random_cell_index(c::VoronoiCellsA, cutoff)
    ntries::Int = 0
    local idx::VoronoiCellIndex
    while true
        pt = random_point(cutoff)
        ntries += 1
        idx = findindex0(c,pt)
        if isvalid(idx) break end
        if ntries > 1000 error("VoronoiCells: can't find a random cell") end
    end
    idx
end


"""
isvalid(idx::VoronoiCellIndex)

Return `true` if `idx` is valid cell index.

## Example

idx = findindex(cells,x,y)
isvalid(idx)

return true if a cell was found containing `(x,y)`, otherwise `false`.
"""
Base.isvalid(iv::VoronoiCellIndex) = iv._ind != 0
splat(iv::VoronoiCellIndex) = (iv._ix,iv._iy,iv._ind)
==(iv1::VoronoiCellIndex, iv2::VoronoiCellIndex) = splat(iv1) == splat(iv2)

"""
    getinvalidcellindex()

return an invalid cell index. This is useful for initializing an index
that will later be checked for validity.
"""
getinvalidcellindex() = VoronoiCellIndex(0,0,0)

getareascale(gcells::VoronoiCellsA) = gcells._areascale
getscale(gcells::VoronoiCellsA) = gcells._scale
getshift(gcells::VoronoiCellsA) = gcells._shift
getcells(gcells::VoronoiCellsA) = gcells._cells
# Note that c[i], c[i,j], c[i,j,k] return very different things
getindex(c::VoronoiCellsA, i::Int) = c._cells[i]
getindex(c::VoronoiCellsA, i::Int, j::Int) = c._grid[i,j]
getindex(c::VoronoiCellsA, i::Int, j::Int, k::Int) =  c._cells[c._grid[i,j][k]]
getindex(c::VoronoiCellsA, iv::VoronoiCellIndex) = getindex(c,splat(iv)...)
endof(c::VoronoiCellsA) = endof(c._cells)
# return integer index into big 1d array of all cells
getcellindex(c::VoronoiCellsA, i::Int, j::Int, k::Int) =  c._grid[i,j][k]
getcellindex(c::VoronoiCellsA, iv::VoronoiCellIndex) = getcellindex(c, splat(iv)...)
Base.length(c::VoronoiCellsA) = length(c._cells)
ngrid(c::VoronoiCellsA) = c._ngrid

# Some of these look broken to me... iscale must take two args
# inverse of scaling of coordinates
iscale(gc::VoronoiCellsA, x) = (x / gc._scale) + gc._shift
iscale(gc::VoronoiCellsA, p::Point2D) = Point2D(iscale(gc,getx(p)),iscale(gc,gety(p)))

Base.scale(gc::VoronoiCellsA, p::Point2D) = Point2D(scale(gc,getx(p)),scale(gc,gety(p)))
Base.scale(gc::VoronoiCellsA, x) = (x - gc._shift) * gc._scale

# for initializing member _grid above
function make_grid_array(ngrid::Int)
    gridcells = Array(Array{Int,1},ngrid,ngrid)
    for i in 1:ngrid
        for j in 1:ngrid
        @inbounds gridcells[i,j] = Array(Int,0)
        end
    end
    gridcells
end

# Find the element in the 2d grid containing the point (x,y)
# Assume input point satisfies 1.0 <= x,y <= 2.0 (actually a bit different)
# Divide this region into ngrid x ngrid grid and return index.
# of box that gp is in.
# This probably throws away carefully preserved precision.
# Should we use width of area rather than 1.0 ?
function find_grid_element(pt::Point2D, ngrid::Int)
    ix::Int = round(Int64, (getx(pt) - 1.0) * ngrid) + 1
    iy::Int = round(Int64, (gety(pt) - 1.0) * ngrid) + 1
    ix > ngrid ? ix = ngrid : nothing
    iy > ngrid ? iy = ngrid : nothing
    ix < 1 ? ix = 1 : nothing
    iy < 1 ? iy = 1 : nothing
    return (ix,iy)  # return indices
end

#find_grid_element(p::Point2D, ngrid::Int) = find_grid_element(getx(p),gety(p),ngrid)
find_grid_element(cell::VoronoiCell, ngrid::Int) = find_grid_element(cell._generator, ngrid)

function make_grid_cells(cells::Array{VoronoiCell,1}, ngrid::Int)
    gridcells = make_grid_array(ngrid)
    n::Int = length(cells)
    @inbounds for i in 1:n
        cell = cells[i]
        (ix::Int,iy::Int) = find_grid_element(cell,ngrid)
        push!(gridcells[ix,iy], i)
    end
    gridcells
end

# BUG!! area scale is wrong here, because cells have been thrown away
# Assign indices of cells in big linear array to arrays in each
# grid element and return new object.
# Cells are assigned to the grid element the generator of the cell is in.
function cellstogrid(cells::Array{VoronoiCell,1}, ngrid::Int)
    gridcells = make_grid_cells(cells,ngrid)
    n = length(cells)
    gshift::Float64 = 1.5
    gscale::Float64 = sqrt(n)
    gareascale::Float64 = convert(Float64,n)
    VoronoiCellsA(ngrid,cells,gridcells, gscale, gshift, gareascale)
end

function cellstogrid(cells::Array{VoronoiCell,1}, ngrid::Int, gscale, gshift, gareascale)
    gridcells = make_grid_cells(cells,ngrid)
    VoronoiCellsA(ngrid,cells,gridcells, gscale, gshift, gareascale)
end

# Version 0.5.0-dev+3385 (2016-04-02 23:53 UTC) throws an error when
# trying to use 'collect' as in commented out code above. So
# we make an explicit loop to create the cell files.
function voronoicells(t::DelaunayTessellation2D, ndiv)
    cells = avoronoicellsnogrid(t)
    ngrid = round(Int,sqrt(length(cells))/ndiv)
    cellstogrid(cells, ngrid)
end

function voronoicells(t::DelaunayTessellation2D)
     voronoicells(t,10) # ndiv defaults to 10
end


# find index in short array of indices into big array of indices of cells
# such that the cell contains point (x,y). If no such index exists in the
# short array, return zero.
# cells -- array of all cells
# indarray -- array of indices into cells corresponding one grid element
# x,y -- point in cell we are searching for
function findindexA(cells::Array{VoronoiCell,1}, indarray::Array{Int,1}, p::Point2D)
    ifound = 0
    @inbounds for i in 1:length(indarray)
        if invoronoicell(cells[indarray[i]],p)
            ifound = i
            break
        end
    end
    ifound
end

# same as findindexA, but gc contains all structures and we identify index (i,j)
# of grid element corresponding to short array of indices.
# Find cell containing (x,y) in bin (i,j) and return index into gc._cells
findindexingrid(gc::VoronoiCellsA, i, j, pt) = findindexA(gc._cells,gc._grid[i,j],pt)

function findindex00(gridcells::VoronoiCellsA, pt::Point2D)
    (ix::Int,iy::Int) = find_grid_element(pt,size(gridcells._grid,1))
    ind::Int = findindexA(gridcells._cells, gridcells._grid[ix,iy], pt)
    return(ix,iy,ind)
end

macro maybe_return_from_find_index()
    esc(quote
        ind = findindexingrid(grc,ix0,iy0,p)
        ind != 0 && return VoronoiCellIndex(ix0,iy0,ind)
        end)
end

# Find index of cell in gc containing point p.
# using the grid makes the search efficient.
# Rounding errors cause point to not be found in a computed grid bin for about 2/3 percent
# of randomly chosen points. In these cases, we look for the point in the neighboring bins.
# Test shows that this works for all random points (no misses found in 10^7 or more trials)
# Each case caught below occurs in 10^6 trials.
function findindex0(grc::VoronoiCellsA, p::Point2D)
    (ix::Int,iy::Int, ind::Int) = findindex00(grc,p)
    ind != 0 && return VoronoiCellIndex(ix,iy,ind)
    if ix > 1
        ix0 = ix-1
        iy0 = iy
        @maybe_return_from_find_index
        if iy > 1
            iy0 = iy-1
            @maybe_return_from_find_index
        end
        if iy < ngrid(grc)
            iy0 = iy+1
            @maybe_return_from_find_index
        end
    end
    if ix < ngrid(grc)
        ix0 = ix+1
        iy0 = iy
        @maybe_return_from_find_index
        if iy > 1
            iy0 = iy-1
            @maybe_return_from_find_index
        end
        if iy < ngrid(grc)
            iy0 = iy+1
            @maybe_return_from_find_index
        end
    end
    if iy > 1
        iy0 = iy-1
        ix0 = ix
        @maybe_return_from_find_index
    end
    if iy < ngrid(grc)
        iy0 = iy+1
        ix0 = ix
        @maybe_return_from_find_index
    end
    return VoronoiCellIndex(ix,iy,ind)
end

# User gives last returned index as hint if new p is close to p for last call
# i.e. assume we are in the same cell. We also assume that the same indices
# ix,iy will be computed. If hint is wrong, we do the usual search.
function findindex0(gridcells::VoronoiCellsA, hint::Int, p::Point2D)
    (ix::Int,iy::Int) = find_grid_element(p,size(gridcells._grid,1))
    length(gridcells[ix,iy]) >= hint && invoronoicell(gridcells[ix,iy,hint],p) && return VoronoiCellIndex(ix,iy,hint)
    findindex0(gridcells,p)
end

findindex(gridcells::VoronoiCellsA, x,y) = findindex0(gridcells,Point2D(x,y))
findindex(gridcells::VoronoiCellsA, p::Point2D) = findindex0(gridcells,p)

findindex(gridcells::VoronoiCellsA, hint::VoronoiCellIndex, p::Point2D) =
    findindex0(gridcells,hint._ind, p)
findindex(gridcells::VoronoiCellsA, hint::VoronoiCellIndex, x,y) = findindex0(gridcells,hint._ind, Point2D(x,y))

# Return false if there is no complete cell containing p
function isexternal(gridcells::VoronoiCellsA, p::Point2D)
    idx = findindex0(gridcells,p)
    idx._ind == 0 && return true
    return false
end

isexternal(gridcells::VoronoiCellsA, x, y) = isexternal(gridcells,Point2D(x,y))

# Careful with these, they don't check for ind == 0
function locate(gridcells::VoronoiCellsA, p::Point2D)
    (ix,iy,ind) = findindex(gridcells,p)
    cellind = gridcells._grid[ix,iy][ind]
    gridcells._cells[cellind]
end

locate(gridcells::VoronoiCellsA,x,y) =  locate(gridcells, Point2D(x,y))

# Return just the tesselation structure.

# Generate a sample of Poisson point process in the plane and return Delaunay tesselation.
# If Distributions is installed, sample the number of generator points npts in the unit square, assuming that
# the mean number is n. Otherwise, we take npts to be n. In any case approxpoissontesselation
# always does the latter.

# We should generate all these with a macro, of course!
if haveDistributions

"""
    poissontesselation(n::Int)

return a tesselation `t::DelaunayTessellation` of a sample of the
Poisson point process in the region of the plane used by `GeometricalPredicates`.
This depends on `Distributions.jl`.
"""    
@eval function poissontesselation{T<:Real}(n::T)
    width = max_coord - min_coord
    npts = rand(Poisson(n))
    a = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:npts]  # this is very fast
    tess = DelaunayTessellation()   # this is 5 or 6 times faster than making the cell structure
    push!(tess,a)
    tess
end
else
 @eval   function poissontesselation(n::Int)
    width = max_coord - min_coord
    a = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:n]  # this is very fast
    tess = DelaunayTessellation()   # this is 5 or 6 times faster than making the cell structure
    push!(tess,a)
    tess
   end
end

# Not an exact sample of the Poisson point process in the plane
"""
    approxpoissontesselation(n::Int)

return a tesselation `t::DelaunayTessellation` of an approximation of a sample of the
Poisson point process in the region of the plane used by `GeometricalPredicates`.
Exactly `n` generator points will be in the sample.
"""
function approxpoissontesselation(n::Int)
    width = max_coord - min_coord
    a = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:n]  # this is very fast
    tess = DelaunayTessellation()   # this is 5 or 6 times faster than making the cell structure
    push!(tess,a)
    tess
end


# Create tesselation of poisson point process sample.
# Return only the efficient cell structure.
# ndiv^2 is the average number of generator points per square in the grid
function poissonvoronoicells(n, ndiv)
    tess = poissontesselation(n)
    gcells = voronoicells(tess,ndiv)
end

# Same as above, but calculate standard grid sized from number of cells
function poissonvoronoicells(n)
    tess = poissontesselation(n)
    gcells = voronoicells(tess)
end

# Generate and return single array of all cells. Do not create grid
# structure. Cannot search efficiently for random points in these cells.
# See below about not using 'collect'
poissonvoronoicellsnogrid(n::Int) = avoronoicellsnogrid(poissontesselation(n))

function approxpoissonvoronoicells(n, ndiv)
    tess = approxpoissontesselation(n)
    gcells = voronoicells(tess,ndiv)
end

function approxpoissonvoronoicells(n)
    tess = approxpoissontesselation(n)
    gcells = voronoicells(tess)
end

approxpoissonvoronoicellsnogrid(n::Int) = avoronoicellsnogrid(approxpoissontesselation(n))

####

# Scaled versions of some functions. The tesselated region, approximately
# 1.0 <= x,y, 2.0  is scaled and shifted to coordinates convenient for the
# user.
slocate(gcells::VoronoiCellsA,x,y) = locate(gcells,iscale(gcells,x),iscale(gcells,y))
slocate(gcells::VoronoiCellsA, p::Point2D) = locate(gridcells, iscale(gcells,p))
sfindindex0(gcells::VoronoiCellsA, x,y) = findindex0(gcells,iscale(gcells,x),iscale(gcells,y))
sfindindex0(gcells::VoronoiCellsA, p::Point2D) = findindex0(gcells, iscale(gcells,p))
sfindindex0(gcells::VoronoiCellsA, hint, x,y) = findindex0(gcells,hint,iscale(gcells,x),iscale(gcells,y))
sfindindex0(gcells::VoronoiCellsA, hint, p::Point2D) = findindex0(gcells, hint, iscale(gcells,p))
sfindindex(gcells::VoronoiCellsA, x,y) = findindex(gcells,iscale(gcells,x),iscale(gcells,y))
sfindindex(gcells::VoronoiCellsA, p::Point2D) = findindex(gcells, iscale(gcells,p))
sfindindex(gridcells::VoronoiCellsA, hint::VoronoiCellIndex, x,y) = findindex(gridcells, hint, iscale(gridcells,x),iscale(gridcells,y))

sarea(gcells::VoronoiCellsA, i::Int) = area(gcells[i]) * getareascale(gcells)
sarea(gcells::VoronoiCellsA, i::Int, j::Int, k::Int) = area(gcells[i,j,k]) * getareascale(gcells)
sarea(gcells::VoronoiCellsA, idx::VoronoiCellIndex) = sarea(gcells, splat(idx)...)

area(gcells::VoronoiCellsA, i::Int, j::Int, k::Int) = area(gcells[i,j,k])
area(gcells::VoronoiCellsA, idx::VoronoiCellIndex) = area(gcells, splat(idx)...)

smaxcoord(gcells) = scale(gcells,max_coord)
smincoord(gcells) = scale(gcells,min_coord)
sgetgenerator(gcells::VoronoiCellsA, c::VoronoiCell) = scale(gcells,c._generator)

#####

function sizeof(c::VoronoiCell)
    s1::Int = sizeof(c._generator)
    s2::Int = sizeof(c._verts)
    s3::Int = 8  # for address of c._verts
    s1 + s2 + s3
end

function sizeof(cells::VoronoiCellsA)
    s1::Int = sizeof(cells._cells)
    s2::Int = sizeof(cells._grid)
    s3::Int = 32 # for four simple members
    s1 + s2 + s3
end

##### Following not exported

# This should be about ndiv^2.
function mean_num_points_per_grid_square(cells::VoronoiCellsA)
    sum::Int = 0
    @inbounds for i in 1:length(cells._grid)
        sum += length(cells._grid[i])
    end
    sum / length(cells._grid)
end

function number_of_generators_in_grid(cells::VoronoiCellsA)
    sum::Int = 0
    @inbounds for i in 1:length(cells._grid)
        sum += length(cells._grid[i])
    end
    sum
end
