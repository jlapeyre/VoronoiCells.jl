using VoronoiDelaunay
using VoronoiCells

# Much of this may not work anymore

# I think getx,gety not implemented for abstract Points
print_pt(str::IOStream, p::GeometricalPredicates.Point2D) = println(str,getx(p), " ", gety(p))

# For plotting triangle, print 4 points
# TODO: check if we can use abstract triangles
function print_triangle(str::IOStream, t::VoronoiDelaunay.DelaunayTriangle)
    print_pt(str,geta(t))
    print_pt(str,getb(t))
    print_pt(str,getc(t))
    print_pt(str,geta(t))
end

function print_delaunay_edges(str::IOStream, t::DelaunayTessellation2D)
    for edge in delaunayedges(t)
        print_pt(str,geta(edge))
        print_pt(str,getb(edge))
        println(str)
    end
end

function print_delaunay_edges(fname::String, t::DelaunayTessellation2D)
    str = open(fname,"w")
    print_delaunay_edges(str, t)
    close(str)    
end

# Print all edges in voronoi cell, one point per line, for plotting
function print_cell(str::IOStream, cell::VoronoiCell)
    verts = cell._verts
    for i in 1:length(verts)
        print_pt(str,verts[i])
    end
    print_pt(str,verts[1])
end

# Print all cells
function print_cells(str::IOStream, cell_array)
    i = 0
    for cell in cell_array
        println(str,"# $i")
        print_cell(str,cell)
        println(str)
        i += 1
    end
end

function print_cells(fname::String,cell_array)
    str = open(fname,"w")
    print_cells(str, cell_array)
    close(str)
end


###### These are for testing code

# prints pairs of common vertices in two triangles
#  a : x
#  b : y
#  c : z
# Where x,y,z are one of a,b,c or nothing.
function print_common_points(tr1::VoronoiDelaunay.DelaunayTriangle,tr2::VoronoiDelaunay.DelaunayTriangle)
    print(" a : ")
    haspt = print_match_point(tr2,geta(tr1))
    if ! haspt println(" x ") end
    print(" b : ")
    haspt = print_match_point(tr2,getb(tr1))
    if ! haspt println(" x ") end
    print(" c : ")
    haspt = print_match_point(tr2,getc(tr1))
    if ! haspt println(" x ") end
end

function print_match_point(tr::VoronoiDelaunay.DelaunayTriangle, p0::Point2D)
    eps0 = eps(min_coord) # same as for max_coord
    has_vertex_f = false    
    pa = geta(tr)
    if distanceL1(pa,p0) <= eps0
        println(" a ")
        has_vertex_f = true
    else
        pb = getb(tr)
        if distanceL1(pb,p0) <= eps0
            println(" b ")
            has_vertex_f = true
        else
            pc = getc(tr)
            if distanceL1(pc,p0) <= eps0
                println(" c ")
                has_vertex_f = true
            end
        end
    end
    has_vertex_f
end
