using VoronoiCells
using Base.Test
using GeometricalPredicates

@test (cells = poissonvoronoicells(10^4) ; true)

cells = poissonvoronoicells(10^2)

@test typeof(sarea(cells,1)) == Float64
@test typeof(area(cells[1])) == Float64

# Test that any point in the bulk is covered by a cell.
function testrand(cells,n)
    c = 0
    s = 1.2
    f = 0.6
    for i in 1:n
        (x,y) = (f*rand()+s,f*rand()+s)
        if VoronoiCells.isexternal(cells, x,y) c += 1 end
    end
    (c,c/n)
end

@test testrand(poissonvoronoicells(10^5),10^5) == (0,0.0)

# Test mean area of cells in the bulk
function meanarea(cells::Array{VoronoiCell,1}, cutoff)
    ma = 0.0
    cnt = 0
    lc = 1 + cutoff
    hc = 2 - cutoff
    for i in 1:length(cells)
        c = cells[i]
        g = getgenerator(c)
        getx(g) < lc && continue
        getx(g) > hc && continue
        gety(g) < lc && continue
        gety(g) > hc && continue
        cnt += 1
        ma += VoronoiCells.area(cells[i])
    end
    ma / cnt
end

@test abs(meanarea(poissonvoronoicells(10^5)._cells,.1) / 1e-5 -1) < 1e-2

