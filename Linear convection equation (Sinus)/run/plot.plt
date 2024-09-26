set terminal pngcairo
set output "result.png"

set xlabel "x"
set ylabel "u"
set yrange [- 1.5 : 1.5]
set grid

plot "result.dat" u 1:2 w l lw 2 lc rgb "red" title "Analytical solution", \
  "" u 1:3 w l lw 2 lc rgb "blue" title "Exact solution", \
  "" u 1:4 w l lw 2 lc rgb "green" title "Numerical solution"
