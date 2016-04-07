#using VoronoiCells
#using GeometricalPredicates

export read_cell_grid, read_n_cell_samples, read_cells, read_one_cellfile

# type of object written to data file
const CELLARRAY = 1
const CELLGRID = 2

type PoissonCellSamples
    _data::Array{Array{VoronoiCell,1},1}
end

Base.length(samples::PoissonCellSamples) = length(samples._data)
Base.getindex(samples::PoissonCellSamples, i) = samples._data[i]

function Base.write(fn::AbstractString, samples::PoissonCellSamples)
    strm = open(fn, "w")
    nsamp = length(samples)
#    write(strm,nsamp)
    for i in 1:nsamp
        cells = samples[1]
        write(strm,cells)
    end
    close(strm)
end

function Base.write(strm::IOStream,cells::Array{VoronoiCell,1})
    ncells = length(cells)
    write(strm,CELLARRAY)
    write(strm,ncells)
    for i in 1:ncells
        cell::VoronoiCell = cells[i]
        gx = getx(cell._generator)
        gy = gety(cell._generator)
        write(strm,gx)
        write(strm,gy)
        nverts = length(cell._verts)
        write(strm,nverts)
        @inbounds for j in 1:nverts
            vert = cell._verts[j]
            write(strm, getx(vert))
            write(strm, gety(vert))
        end
    end
end

function Base.write(strm::IOStream,gcells::VoronoiCellsA)
    write(strm,CELLGRID)
    write(strm,gcells._ngrid)
    write(strm,gcells._scale)
    write(strm,gcells._shift)
    write(strm,gcells._areascale)
    write(strm, gcells._cells)
end

function read_cell_grid(strm::IOStream)
    buf = zeros(Int,1)
    dbuf = zeros(Float64,3)
    read!(strm,buf)
    ngrid = buf[1]
    read!(strm,dbuf)
    scale = dbuf[1]
    shift = dbuf[2]
    areascale = dbuf[3]
    read!(strm,buf)
    data_type = buf[1]
    data_type == CELLARRAY || error("Unknown data type $data_type when reading cell grid")
    cells = read_cells(strm)
    gcells = cellstogrid(cells,ngrid, scale, shift, areascale)
    gcells
end


function Base.write(fn::AbstractString, cells::Array{VoronoiCell,1})
    samples = PoissonCellSamples(Array(Array{VoronoiCell,1},0))
    push!(samples._data,cells)
#    write_cells(fn,samples)
    write(fn,samples)
end

function Base.write(fn::AbstractString, gcells::VoronoiCellsA)
    strm = open(fn, "w")
    write(strm,gcells)
    close(strm)
end

function read_n_cell_samples(strm::IOStream)
    buf = zeros(Int,1)
    read!(strm,buf)
    nsamp = buf[1]
end

function read_cells(strm::IOStream)
    buf = zeros(Int,1)
    read!(strm,buf)
    ncells = buf[1]
    cells = Array(VoronoiCell,0)
    dbuf = zeros(Float64,1)
    for i in 1:ncells
        read!(strm,dbuf)
        gx = dbuf[1]
        read!(strm,dbuf)
        gy = dbuf[1]
        verts = Array(Point2D,0)
        read!(strm,buf)
        nverts = buf[1]
        for j in 1:nverts
            read!(strm,dbuf)
            vx = dbuf[1]
            read!(strm,dbuf)
            vy = dbuf[1]
            push!(verts, Point2D(vx,vy))
        end
        push!(cells,VoronoiCell(Point2D(gx,gy),verts))
    end
    cells
end

function read_cells_top(strm::IOStream)
    buf = zeros(Int,1)
    read!(strm,buf)
    data_type = buf[1]
    if data_type == CELLARRAY
        return read_cells(strm)
    elseif data_type == CELLGRID
        return read_cell_grid(strm)
    else
        error("read_cells_top: Unrecognized data type ", data_type)
    end
end

function read_one_cellfile(fn::AbstractString)
    strm = open(fn,"r")
#    nsamps = read_n_cell_samples(strm)
#    println("n samples ", nsamps)
#    read_cells(strm)
    read_cells_top(strm)
end

function ==(c1::VoronoiCell, c2::VoronoiCell)
    c1._generator == c2._generator || return false
    n1 = length(c1._verts)
    n2 = length(c1._verts)
    n1 == n2 || return false
    for i in 1:n1
        c1._verts[i] == c2._verts[i] || return false
    end
    true
end

# This can be used to check for equality of written and read cells with grid

function ==(cells1::VoronoiCellsA, cells2::VoronoiCellsA)
    length(cells1) != length(cells2) && return false
    for i in 1:length(cells1)
        cells1[i] != cells2[i] && return false
    end
    getshift(cells1) != getshift(cells2) && return false
    getscale(cells1) != getscale(cells2) && return false
    getareascale(cells1) != getareascale(cells2) && return false
    true
end
