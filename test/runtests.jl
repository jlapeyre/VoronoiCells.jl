using VoronoiCells
using Base.Test

@test (cells = poissonvoronoicells(10^4) ; true)

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
