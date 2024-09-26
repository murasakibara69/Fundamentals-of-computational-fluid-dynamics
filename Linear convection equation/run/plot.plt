set terminal png
set output "result.png"

set grid
set xlabel "x"
set ylabel "u"
set xrange [0 : 1]
set yrange [0 : 0.8]

plot "result.dat" u 1:2 w l lw 2 lc rgb "red" title "Exact solution", "" u 1:3 w l lw 2 lc rgb "blue" title "Numerical solution"
