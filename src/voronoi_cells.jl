#module VoronoiCells

using VoronoiDelaunay
using GeometricalPredicates
import GeometricalPredicates
import GeometricalPredicates: area  # needed so that area() can be called in runtests.jl
import VoronoiDelaunay.isexternal, VoronoiDelaunay.locate, VoronoiDelaunay.findindex
import Base: getindex, ==, -

export VoronoiCell, VoronoiCellsA, voronoicells, voronoicellsnogrid, voronoicells2, findindex, locate, invoronoicell, area, avoronoicellsnogrid
export getcellindex,isexternal, nverts, nedges, getgenerator, scale, iscale

export VoronoiCellIdx, isvalid, getinvalidcellindex

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
# P(x,y) lies
lineside(x1,y1,x2,y2,x,y) = (y-y1)*(x2-x1) - (x-x1)*(y2-y1)

# true if p is inside poly, ie. if lineside gives the same sign for all edges in Voronoi cell
function inconvexpolygon(poly::Array{Point2D,1}, x,y)
    is_inconvexpolygon = true
    n = length(poly)
    p1 = poly[n]
    p2 = poly[1]
    s = lineside(getx(p1),gety(p1),getx(p2),gety(p2),x,y)
    for i in 1:length(poly)-1
        p1 = poly[i]
        p2 = poly[i+1]
        s1 = lineside(getx(p1),gety(p1),getx(p2),gety(p2),x,y)
        if sign(s1) != sign(s)
            is_inconvexpolygon = false
            break
        end
    end
    is_inconvexpolygon
end
inconvexpolygon(poly::Array{Point2D,1}, p::Point2D) = inconvexpolygon(poly,getx(p),gety(p))

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
function get_neighbor(trigs,trig,iv)
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
function get_neighbor_index(trigs,trig,iv)
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
    elseif v1 == 3
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
function next_cell_triangle(trigs,tr,iv_gen::Int,iv_opp::Int, iv_other::Int)
    eps0 = eps(min_coord)
    tr2 = get_neighbor(trigs,tr,iv_opp)
    itr2 = get_neighbor_index(trigs,tr,iv_opp)
    iv_gen2 = match_vertex(tr2,get_vertex(tr,iv_gen))
    iv_opp2 = match_vertex(tr2,get_vertex(tr,iv_other))
    iv_other2 = match_other_vertex(iv_gen2, iv_opp2)
    return (tr2, itr2, iv_gen2, iv_opp2, iv_other2)
end

#### Begin code specific to VoronoiCells

# 2D
# _generator is one of the points in the point process. There is exactly one such point in each
# cell.
type VoronoiCell
    _generator::Point2D
    _verts::Array{Point2D,1}
    # neigbhors ?
end

VoronoiCell(generator) = VoronoiCell(generator, Array(Point2D, 0))
nverts(c::VoronoiCell) = length(c._verts)
nedges(c::VoronoiCell) = length(c._verts) + 1
getgenerator(c::VoronoiCell) = c._generator

# Area of irregular polygon. Don't make use of convex property.
function area(c::VoronoiCell)
    vs = c._verts
    A = 0.0
    n = length(vs)
    for i in 1:n-1
        A += getx(vs[i])*gety(vs[i+1]) -  getx(vs[i+1])*gety(vs[i])
    end
    A += getx(vs[n])*gety(vs[1]) -  getx(vs[1])*gety(vs[n])
    0.5 * abs(A)
end

# is any vertex on boundary of cell outside of allowed range
function isexternal(c::VoronoiCell)
    vs = c._verts
    n = length(vs)
    found_external = false
    for i in 1:n
        xc = getx(vs[i])
        if xc < min_coord || xc > max_coord
            found_external = true
            break
        end
        yc = gety(vs[i])
        if yc < min_coord || yc > max_coord
            found_external = true
            break
        end
    end
    found_external
end

invoronoicell(c::VoronoiCell, p::Point2D) = inconvexpolygon(c._verts,p)
invoronoicell(c::VoronoiCell, x,y) = inconvexpolygon(c._verts,x,y)

# Find the index in an array of cells of the (first) cell containing point p,
# or zero if no cell contains p.
# Use linear search. We do not use this function on type VoronoiCellsA below.
# It is too expensive to search all cells. Searching for one point in each of 10^6
# cells takes on the order of 1 day cpu time.
# Below, there is a method for this function for VoronoiCellsA.
function findindex(cells::Array{VoronoiCell,1}, p::Point2D)
    ifound = 0
    for i in 1:length(cells)
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
#        if (visited[itr2] || isexternal(tr2)) && tr2 != tr0
        if visited[itr2] && tr2 != tr0
            is_visited = true
            break
        end
        # if isexternal(tr)
        #     is_visited = true
        #     break
        # end
        push!(cell._verts, circumcenter(tr))
        if tr2 == tr0 break end
        tr = tr2
    end
    (is_visited,cell)
end

# function voronoicells2(t::DelaunayTessellation2D)
#     visited = zeros(Bool, t._last_trig_index)
#     visited[1] = true
#     function voronoicelliterator()
#         j = 0
# 	for ix in 2:t._last_trig_index
# 	    visited[ix] && continue
# 	    const tr = t._trigs[ix]
# 	    visited[ix] = true
#             isexternal(tr) && continue
#             for iv_gen in 1:3

#                 iv_opp = mod1(iv_gen+1,3) # pick one of the other two vertices

#                 # ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
#                 # visited[ix2] && continue

#                 #                isexternal(t._trigs[ix2]) && continue

#                 # iv_opp = mod1(iv_gen+2,3) # pick the remaining vertex
#                 # ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
#                 # visited[ix2] && continue

#                 #                isexternal(t._trigs[ix2]) && continue

#                 (is_visited, cell) = find_cell(t._trigs, tr, iv_gen, iv_opp, visited)
#                 is_visited && continue
#                 j += 1
# #                j > 10000 && return
#                 isexternal(cell) && continue # just checking each triangle does not seem to work.
#                 produce(cell)
#             end
#         end
#     end
#     Task(voronoicelliterator)
# end

function voronoicellsnogrid(t::DelaunayTessellation2D)
    visited = zeros(Bool, t._last_trig_index)
    visited[1] = true
    function voronoicelliterator()
        j = 0
	for ix in 2:t._last_trig_index
	    visited[ix] && continue
	    const tr = t._trigs[ix]
	    visited[ix] = true
            isexternal(tr) && continue
            for iv_gen in 1:3

                iv_opp = mod1(iv_gen+1,3) # pick one of the other two vertices
                ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
                visited[ix2] && continue
                #                isexternal(t._trigs[ix2]) && continue

                iv_opp = mod1(iv_gen+2,3) # pick the remaining vertex
                ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
                visited[ix2] && continue

                #                isexternal(t._trigs[ix2]) && continue

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

function avoronoicellsnogrid(t::DelaunayTessellation2D)
    visited = zeros(Bool, t._last_trig_index)
    visited[1] = true
    cells = Array(VoronoiCell,0)
    j = 0
    for ix in 2:t._last_trig_index
	visited[ix] && continue
	const tr = t._trigs[ix]
	visited[ix] = true
        isexternal(tr) && continue
        for iv_gen in 1:3

            iv_opp = mod1(iv_gen+1,3) # pick one of the other two vertices
            ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
            visited[ix2] && continue
            #                isexternal(t._trigs[ix2]) && continue

            iv_opp = mod1(iv_gen+2,3) # pick the remaining vertex
            ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
            visited[ix2] && continue

            #                isexternal(t._trigs[ix2]) && continue

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

export cellstogrid, poissontesselation, poissonvoronoicells, poissonvoronoicellsnogrid,
       ngrid, findindex0
export slocate, sfindindex0, sfindindex, sarea, smaxcoord, smincoord, sgetgenerator
export getareascale, getscale, getshift

#### VoronoiCellsA

# Voronoi cells arranged in efficiently searchable array
# structure is _ngrid by _ngrid array of arrays of indices into the linear array.
# elements in the 2d array correspond to square regions covering area of point process.
#  of all VoronoiCell's.
# _cells -- linear array of all cells (currently in (almost) no order)
# _grid -- 2d array of arrays of indices into _cells.
# _scale, _shift, _areascale. For convenience of user. See scaled versions of functions below.
type VoronoiCellsA
    _ngrid::Int
    _cells::Array{VoronoiCell,1}
    _grid::Array{Array{Int,1},2}
    _scale::Float64
    _shift::Float64
    _areascale::Float64
end

# index into grid and array structure VoronoiCellsA
immutable VoronoiCellIdx
    _ix::Int
    _iy::Int
    _ind::Int
end


Base.isvalid(iv::VoronoiCellIdx) = iv._ind != 0
splat(iv::VoronoiCellIdx) = (iv._ix,iv._iy,iv._ind)
==(iv1::VoronoiCellIdx, iv2::VoronoiCellIdx) = splat(iv1) == splat(iv2)
# ==(iv1::VoronoiCellIdx, iv2::VoronoiCellIdx) = iv1._ix == iv2._ix && iv1._iy == iv2._iy && iv1._ind == iv2._ind
getinvalidcellindex() = VoronoiCellIdx(0,0,0)

getareascale(gcells::VoronoiCellsA) = gcells._areascale
getscale(gcells::VoronoiCellsA) = gcells._scale
getshift(gcells::VoronoiCellsA) = gcells._shift
# Note that c[i], c[i,j], c[i,j,k] return very different things
getindex(c::VoronoiCellsA, i::Int) = c._cells[i]
getindex(c::VoronoiCellsA, i::Int, j::Int) = c._grid[i,j]
getindex(c::VoronoiCellsA, i::Int, j::Int, k::Int) =  c._cells[c._grid[i,j][k]]
getindex(c::VoronoiCellsA, iv::VoronoiCellIdx) = getindex(c,splat(iv)...)
#getindex(c::VoronoiCellsA, iv::VoronoiCellIdx) = getindex(c,iv._ix, iv._iy, iv._ind)
# return integer index into big 1d array of all cells
getcellindex(c::VoronoiCellsA, i::Int, j::Int, k::Int) =  c._grid[i,j][k]
getcellindex(c::VoronoiCellsA, iv::VoronoiCellIdx) = getcellindex(c, splat(iv)...)
#getcellindex(c::VoronoiCellsA, iv::VoronoiCellIdx) = getcellindex(c, iv._ix, iv._iy, iv._ind)
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
            gridcells[i,j] = Array(Int,0)
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
function find_grid_element(x::Float64, y::Float64, ngrid::Int)
    ix = round(Int64, (x - 1.0) * ngrid) + 1
    iy = round(Int64, (y - 1.0) * ngrid) + 1
    ix > ngrid ? ix = ngrid : nothing
    iy > ngrid ? iy = ngrid : nothing
    ix < 1 ? ix = 1 : nothing
    iy < 1 ? iy = 1 : nothing
    return (ix,iy)  # return indices
end

find_grid_element(p::Point2D, ngrid::Int) = find_grid_element(getx(p),gety(p),ngrid)
find_grid_element(cell::VoronoiCell, ngrid::Int) = find_grid_element(cell._generator, ngrid)

# Assign indices of cells in big linear array to arrays in each
# grid element and return new object.
# Cells are assigned to the grid element the generator of the cell is in.
function cellstogrid(cells::Array{VoronoiCell,1}, ngrid::Int)
    gridcells = make_grid_array(ngrid)
    for i in 1:length(cells)
        cell = cells[i]
        (ix,iy) = find_grid_element(cell,ngrid)
        push!(gridcells[ix,iy], i)
    end
    VoronoiCellsA(ngrid,cells,gridcells,1.0,0.0,1.0)
end

# generate big array of cells from tesselation and store them in grid structure
function voronoicells(t::DelaunayTessellation2D, ngrid::Int)
    cells = collect(VoronoiCell,voronoicellsnogrid(t))
    cellstogrid(cells, ngrid)
end


# generate as above, but choose a good default number of elements in the grid based on the
# number of cells (good if they are uniformly distributed)
# function voronoicells(t::DelaunayTessellation2D)
#     cells = collect(VoronoiCell,voronoicellsnogrid(t))
#     ngrid = round(Int,sqrt(length(cells))/10)
#     cellstogrid(cells, ngrid)
# end

# Version 0.5.0-dev+3385 (2016-04-02 23:53 UTC) throws an error when
# trying to use 'collect' as in commented out code above. So
# we make an explicit loop to create the cell files.
function voronoicells(t::DelaunayTessellation2D)
    celltask = voronoicellsnogrid(t)
    cells = Array(VoronoiCell,0)
    for cell in celltask
        push!(cells, cell)
    end
    ngrid = round(Int,sqrt(length(cells))/10)
    cellstogrid(cells, ngrid)
end

# find index in short array of indices into big array of indices of cells
# such that the cell contains point (x,y). If no such index exists in the
# short array, return zero.
# cells -- array of all cells
# indarray -- array of indices into cells corresponding one grid element
# x,y -- point in cell we are searching for
function findindexA(cells::Array{VoronoiCell,1}, indarray::Array{Int,1}, x,y)
    ifound = 0
    for i in 1:length(indarray)
        if invoronoicell(cells[indarray[i]],x,y)
            ifound = i
            break
        end
    end
    ifound
end

# Find index into gc._cells of cell containing (x,y) searching in indarray
findindexA(cells::Array{VoronoiCell,1}, indarray::Array{Int,1}, p::Point2D) = findindexA(cells,indarray,getx(p),gety(p))

# same as findindexA, but gc contains all structures and we identify index (i,j)
# of grid element corresponding to short array of indices.
# Find cell containing (x,y) in bin (i,j) and return index into gc._cells
findindexingrid(gc::VoronoiCellsA, i, j, x, y) = findindexA(gc._cells,gc._grid[i,j],x,y)


# TODO: make cleaner use of (x,y)  <--> Point2D(x,y)
function findindex00(gridcells::VoronoiCellsA, x::Float64, y::Float64)
    (ix::Int,iy::Int) = find_grid_element(x,y,size(gridcells._grid,1))
    ind::Int = findindexA(gridcells._cells, gridcells._grid[ix,iy], Point2D(x,y))
    return(ix,iy,ind)
end

# Find index of cell in gc containing point p.
# using the grid makes the search efficient.
# Rounding errors cause point to not be found in a computed grid bin for about 2/3 percent
# of randomly chosen points. In these cases, we look for the point in the neighboring bins.
# Test shows that this works for all random points (no misses found in 10^7 or more trials)
# Each case below occurs in 10^6 trials.
#function findindex0(grc::VoronoiCellsA, p::Point2D)
function findindex0(grc::VoronoiCellsA, x, y)
#    (x,y) = (getx(p),gety(p))
    (ix::Int,iy::Int, ind::Int) = findindex00(grc,x,y)
    ind != 0 && return (ix,iy,ind)
    if ix > 1
        ix0 = ix-1
        iy0 = iy
        ind = findindexingrid(grc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
        if iy > 1
            iy0 = iy-1
            ind = findindexingrid(grc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
        if iy < ngrid(grc)
            iy0 = iy+1
            ind = findindexingrid(grc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
    end
    if ix < ngrid(grc)
        ix0 = ix+1
        iy0 = iy
        ind = findindexingrid(grc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
        if iy > 1
            iy0 = iy-1
            ind = findindexingrid(grc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
        if iy < ngrid(grc)
            iy0 = iy+1
            ind = findindexingrid(grc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
    end
    if iy > 1
        iy0 = iy-1
        ix0 = ix
        ind = findindexingrid(grc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
    end
    if iy < ngrid(grc)
        iy0 = iy+1
        ix0 = ix
        ind = findindexingrid(grc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
    end
    return (ix,iy,ind)
end

# User gives last returned index as hint if new p is close to p for last call
# i.e. assume we are in the same cell. We also assume that the same indices
# ix,iy will be computed. If hint is wrong, we do the usual search.
#function findindex0(gridcells::VoronoiCellsA, hint::Int, p::Point2D)
function findindex0(gridcells::VoronoiCellsA, hint::Int, x, y)
    (ix,iy) = find_grid_element(x,y,size(gridcells._grid,1))
    length(gridcells[ix,iy]) >= hint && invoronoicell(gridcells[ix,iy,hint],x,y) && return (ix,iy,hint)
    findindex0(gridcells,x,y)
end
# function findindex0(gridcells::VoronoiCellsA, hint::Int, p::Point2D)
#     (ix,iy) = find_grid_element(p,size(gridcells._grid,1))
#     length(gridcells[ix,iy] >= hint) && invoronoicell(gridcells[ix,iy,hint],p) && return (ix,iy,hint)
#     findindex0(gridcells,p)
# end

#findindex0(gridcells::VoronoiCellsA, x,y) = findindex0(gridcells, Point2D(x,y))
findindex0(gridcells::VoronoiCellsA, p::Point2D) = findindex0(gridcells, getx(p), gety(p))
#findindex0(gridcells::VoronoiCellsA, hint::Int, x,y) = findindex0(gridcells, hint, Point2D(x,y))
findindex0(gridcells::VoronoiCellsA, hint, p::Point2D) = findindex0(gridcells, hint, getx(p), gety(p))
findindex(gridcells::VoronoiCellsA, x,y) = VoronoiCellIdx(findindex0(gridcells,x,y)...)
findindex(gridcells::VoronoiCellsA, hint::VoronoiCellIdx, x,y) = VoronoiCellIdx(findindex0(gridcells,hint._ind, x,y)...)

# Probably don't use this, because we would have to prevent errors first,
# which is expensive.
function findindex1(gridcells::VoronoiCellsA, p::Point2D)
    (ix,iy,ind) = findindex0(gridcells,p)
    ind == 0 && error("locate: Can't find grid box for point ", p, ".")
    return(ix,iy,ind)
end
findindex1(gridcells::VoronoiCellsA, x,y) = findindex(gridcells, Point2D(x,y))

# Return false if there is no complete cell containing p
function isexternal(gridcells::VoronoiCellsA, p::Point2D)
    (ix,iy,ind) = findindex0(gridcells,p)
    ind == 0 && return true
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
function poissontesselation(n::Int)
    width = max_coord - min_coord
    a = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:n]
    tess = DelaunayTessellation()
    push!(tess,a)
    tess
end

# Create tesselation of poisson point process sample.
# Return only the efficient cell structure.
function poissonvoronoicells(n::Int,ngrid::Int)
    tess = poissontesselation(n)
    gcells = voronoicells(tess,ngrid)
    standard_scale_and_shift!(gcells,n)
    gcells
end

# Same as above, but calculate standard grid sized from number of cells
function poissonvoronoicells(n::Int)
    tess = poissontesselation(n)
    gcells = voronoicells(tess)
    standard_scale_and_shift!(gcells,n)
    gcells
end

# Generate and return single array of all cells. Do not create grid
# structure. Cannot search efficiently for random points in these cells.
function poissonvoronoicellsnogrid(n::Int)
    tess = poissontesselation(n)
    cells = collect(VoronoiCell,voronoicellsnogrid(tess))
end

####

# For Poisson point process, origin is at zero zero and average cell
# size is 1.
function standard_scale_and_shift!(gcells::VoronoiCellsA, n::Int)
    gcells._shift = 1.5
    gcells._scale = sqrt(n)
    gcells._areascale = convert(Float64,n)
    nothing
end

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
sfindindex(gridcells::VoronoiCellsA, hint::VoronoiCellIdx, x,y) = findindex(gridcells, hint, iscale(gridcells,x),iscale(gridcells,y))
sarea(gcells::VoronoiCellsA, i::Int) = area(gcells[i]) * getareascale(gcells)
sarea(gcells::VoronoiCellsA, i::Int, j::Int, k::Int) = area(gcells[i,j,k]) * getareascale(gcells)
sarea(gcells::VoronoiCellsA, idx::VoronoiCellIdx) = sarea(gcells, splat(idx)...)
#sarea(gcells::VoronoiCellsA, idx::VoronoiCellIdx) = sarea(gcells, idx._ix, idx._iy, idx._ind)
smaxcoord(gcells) = scale(gcells,max_coord)
smincoord(gcells) = scale(gcells,min_coord)
sgetgenerator(gcells::VoronoiCellsA, c::VoronoiCell) = scale(gcells,c._generator)

#end # module
