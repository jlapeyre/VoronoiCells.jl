module VoronoiCells

using VoronoiDelaunay
using GeometricalPredicates
import GeometricalPredicates
import VoronoiDelaunay.isexternal, VoronoiDelaunay.locate, VoronoiDelaunay.findindex

export VoronoiCell, voronoicells, voronoicellsnogrid, voronoicells2, findindex, locate, invoronoicell, area
export isexternal

const min_coord = GeometricalPredicates.min_coord + eps(Float64)
const max_coord = GeometricalPredicates.max_coord - eps(Float64)

#### Some functions that only operate on GeometricalPredicates

# Euclidean distance
function distance(p1::Point2D,p2::Point2D)
    (getx(p2)-getx(p1))^2 + (gety(p2)-gety(p1))^2
end

# Maybe faster than Euclidean, and is sufficient for use below
function distanceL1(p1::Point2D,p2::Point2D)
    abs(getx(p2)-getx(p1)) + abs(gety(p2)-gety(p1))
end

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

# Note regarding DelaunayTriangle:
# tr._neighbour_a is the neighboring triangle sharing the edge in tr that is opposite point a
# tr._neighbour_b, tr._neighbour_c contain the vertex a, while tr._neighbour_a does not.

#### Some functions that only operate on VoronoiDelaunay
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

#### Begin code specific to VoronoiCells

# 2D
type VoronoiCell
    _generator
    _verts::Array{Point2D,1}
    # neigbhors ?
end

function VoronoiCell(generator)
    VoronoiCell(generator, Array(Point2D, 0))
end

function npoints(c::VoronoiCell)
    length(c._verts)
end

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

# linear search
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
function locate(cells::Array{VoronoiCell,1}, p::Point2D)
    cells[findindex(cells,p)]
end

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

function next_cell_triangle(trigs,tr,iv_gen,iv_opp, iv_other)
    eps0 = eps(min_coord)
    tr2 = get_neighbor(trigs,tr,iv_opp)
    itr2 = get_neighbor_index(trigs,tr,iv_opp)
    iv_gen2 = match_vertex(tr2,get_vertex(tr,iv_gen))
    iv_opp2 = match_vertex(tr2,get_vertex(tr,iv_other))
    iv_other2 = match_other_vertex(iv_gen2, iv_opp2)
    return (tr2, itr2, iv_gen2, iv_opp2, iv_other2)
end

function find_cell(trigs,tr, iv_gen, iv_opp, visited)
    cell = VoronoiCell(get_vertex(tr, iv_gen))
    iv_other = match_other_vertex(iv_gen,iv_opp)
    tr0 = tr
    is_visited = false
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

function voronoicells2(t::DelaunayTessellation2D)
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

                # ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
                # visited[ix2] && continue

                #                isexternal(t._trigs[ix2]) && continue

                # iv_opp = mod1(iv_gen+2,3) # pick the remaining vertex
                # ix2 = get_neighbor_index(t._trigs,tr,iv_opp)
                # visited[ix2] && continue

                #                isexternal(t._trigs[ix2]) && continue

                (is_visited, cell) = find_cell(t._trigs, tr, iv_gen, iv_opp, visited)
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

                (is_visited, cell) = find_cell(t._trigs, tr, iv_gen, iv_opp, visited)
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

export cellstogrid, poissontesselation, poissonvoronoicells, ngrid, findindex0
export slocate, sfindindex0, sfindindex, sarea
export getareascale, getscale, getshift

# Voronoi cells arranged in efficiently searchable array
type VoronoiCellsA
    _ngrid::Int
    _cells::Array{VoronoiCell,1}
    _grid::Array{Array{Int,1},2}
    _scale::Float64
    _shift::Float64
    _areascale::Float64
end

getareascale(gcells::VoronoiCellsA) = gcells._areascale
getscale(gcells::VoronoiCellsA) = gcells._scale
getshift(gcells::VoronoiCellsA) = gcells._shift
getindex(c::VoronoiCellsA, i::Int) = c._cells[i]
getindex(c::VoronoiCellsA, i::Int, j::Int) = c._grid[i,j]
getindex(c::VoronoiCellsA, i::Int, j::Int, k::Int) =  c._cells[c._grid[i,j][k]]
Base.length(c::VoronoiCellsA) = length(c._cells)
ngrid(c::VoronoiCellsA) = c._ngrid

# inverse of scaling of coordinates
iscale(gc::VoronoiCellsA, x) = (x / gc._scale) + gc._shift
iscale(gc::VoronoiCellsA, p::Point2D) = Point2D(iscale(getx(p)),iscale(gety(p)))

function make_grid_array(ngrid::Int)
    gridcells = Array(Array{Int,1},ngrid,ngrid)
    for i in 1:ngrid
        for j in 1:ngrid
            gridcells[i,j] = Array(Int,0)
        end
    end
    gridcells
end

# Assume input point satisfies 1.0 <= x,y <= 2.0 (actually a bit different)
# Divide this region into ngrid x ngrid grid and return index
# of box that gp is in.
# This probably throws away carefully preserved precision
function find_generator_bin(x::Float64, y::Float64, ngrid::Int)
    ix = round(Int64, (x - 1.0) * ngrid) + 1
    iy = round(Int64, (y - 1.0) * ngrid) + 1
    ix > ngrid ? ix = ngrid : nothing
    iy > ngrid ? iy = ngrid : nothing
    ix < 1 ? ix = 1 : nothing
    iy < 1 ? iy = 1 : nothing
    return (ix,iy)
end

find_generator_bin(p::Point2D, ngrid::Int) = find_generator_bin(getx(p),gety(p),ngrid)
find_generator_bin(cell::VoronoiCell, ngrid::Int) = find_generator_bin(cell._generator, ngrid)

function cellstogrid(cells::Array{VoronoiCell,1}, ngrid::Int)
    gridcells = make_grid_array(ngrid)
    for i in 1:length(cells)
        cell = cells[i]
        (ix,iy) = find_generator_bin(cell,ngrid)
        push!(gridcells[ix,iy], i)
    end
    VoronoiCellsA(ngrid,cells,gridcells,1.0,0.0,1.0)
end

function voronoicells(t::DelaunayTessellation2D, ngrid::Int)
    cells = collect(VoronoiCell,voronoicellsnogrid(t))
    cellstogrid(cells, ngrid)
end

function voronoicells(t::DelaunayTessellation2D)
    cells = collect(VoronoiCell,voronoicellsnogrid(t))    
    ngrid = round(Int,sqrt(length(cells))/10)
    cellstogrid(cells, ngrid)
end

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

# Find cell containing (x,y) in bin (i,j) and return index into gc._cells
findindexingrid(gc::VoronoiCellsA, i, j, x, y) = findindexA(gc._cells,gc._grid[i,j],x,y)

#findindexA(cells::Array{VoronoiCell,1}, indarray::Array{Int,1}, p::Point2D) = findindexA(cells,indarray,getx(p),gety(p))

function findindex(gridcells::VoronoiCellsA, p::Point2D)
    (ix,iy,ind) = findindex0(gridcells,p)
    ind == 0 && error("locate: Can't find grid box for point ", p, ".")
    return(ix,iy,ind)
end

findindex(gridcells::VoronoiCellsA, x,y) = findindex(gridcells, Point2D(x,y))

# TODO: make cleaner use of (x,y)  <--> Point2D(x,y)
function findindex00(gridcells::VoronoiCellsA, x::Float64, y::Float64)
    (ix::Int,iy::Int) = find_generator_bin(x,y,size(gridcells._grid,1))
    ind::Int = findindexA(gridcells._cells, gridcells._grid[ix,iy], Point2D(x,y))
    return(ix,iy,ind)
end

# Rounding errors cause point to not be found in a computed grid bin for about 2/3 percent
# of randomly chosen points. In these cases, we look for the point in the neighboring bins.
# Test shows that this works for all random points (no misses found in 10^7 or more trials)
function findindex0(gc::VoronoiCellsA, p::Point2D)
    (x,y) = (getx(p),gety(p))
    (ix::Int,iy::Int, ind::Int) = findindex00(gc,x,y)
    ind != 0 && return (ix,iy,ind)
    if ix > 1
        ix0 = ix-1
        iy0 = iy
        ind = findindexingrid(gc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
        if iy > 1
            iy0 = iy-1
            ind = findindexingrid(gc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
        if iy < ngrid(gc)
            iy0 = iy+1
            ind = findindexingrid(gc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
    end
    if ix < ngrid(gc)
        ix0 = ix+1
        iy0 = iy
        ind = findindexingrid(gc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
        if iy > 1
            iy0 = iy-1
            ind = findindexingrid(gc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
        if iy < ngrid(gc)
            iy0 = iy+1
            ind = findindexingrid(gc,ix0,iy0,x,y)
            ind != 0 && return (ix0,iy0,ind)
        end
    end
    if iy > 1
        iy0 = iy-1
        ix0 = ix
        ind = findindexingrid(gc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
    end
    if iy < ngrid(gc)
        iy0 = iy+1
        ix0 = ix
        ind = findindexingrid(gc,ix0,iy0,x,y)
        ind != 0 && return (ix0,iy0,ind)
    end    

    return (ix,iy,ind)
end

# User gives last returned index as hint if new p is close to p for last call
function findindex0(gridcells::VoronoiCellsA, hint::Int, p::Point2D)
    (ix,iy) = find_generator_bin(p,size(gridcells._grid,1))
    invoronoicell(gridcells[ix,iy,hint],p) && return (ix,iy,hint)
    findindex0(gridcells,p)
end

findindex0(gridcells::VoronoiCellsA, x,y) = findindex0(gridcells, Point2D(x,y))
findindex0(gridcells::VoronoiCellsA, hint::Int, x,y) = findindex0(gridcells, hint, Point2D(x,y))

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

function poissontesselation(n::Int)
    width = max_coord - min_coord
    a = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:n]
    tess = DelaunayTessellation()
    push!(tess,a)
    tess
end

function poissonvoronoicells(n::Int,ngrid::Int)
    tess = poissontesselation(n)
    gcells = voronoicells(tess,ngrid)
    standard_scale_and_shift!(gcells,n)
    gcells    
end

function poissonvoronoicells(n::Int)
    tess = poissontesselation(n)
    gcells = voronoicells(tess)
    standard_scale_and_shift!(gcells,n)
    gcells
end

####

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
slocate(gcells::VoronoiCellsA, p::Point2D) = locate(gridcells, iscale(p))
sfindindex0(gcells::VoronoiCellsA, x,y) = findindex0(gcells,iscale(gcells,x),iscale(gcells,y))
sfindindex0(gcells::VoronoiCellsA, p::Point2D) = findindex0(gcells, iscale(p))
sfindindex0(gcells::VoronoiCellsA, hint, x,y) = findindex0(gcells,hint,iscale(gcells,x),iscale(gcells,y))
sfindindex0(gcells::VoronoiCellsA, hint, p::Point2D) = findindex0(gcells, hint, iscale(p))
sfindindex(gcells::VoronoiCellsA, x,y) = findindex(gcells,iscale(gcells,x),iscale(gcells,y))
sfindindex(gcells::VoronoiCellsA, p::Point2D) = findindex(gcells, iscale(p))
sarea(gcells::VoronoiCellsA, i::Int) = area(gcells[i]) * getareascale(gcells)
sarea(gcells::VoronoiCellsA, i::Int, j::Int, k::Int) = area(gcells[i,j,k]) * getareascale(gcells)

end # module
