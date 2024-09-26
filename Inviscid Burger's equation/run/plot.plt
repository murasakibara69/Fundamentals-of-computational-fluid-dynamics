set terminal pngcairo
set output "Result.png"

set xrange [0 : 1]
set yrange [0 : 0.8]
set xlabel "x"
set ylabel "u"

set grid

plot "result.dat" u 1:2 w l lw 2 lc rgb "red" title "Analytical solution", "" u 1:3 w l lw 2 lc rgb "blue" title "Numerical solution"
