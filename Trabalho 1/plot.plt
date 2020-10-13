plot "output_rk4.dat" using 2:4 title "Y" with lines, "output_rk4.dat" using 2:5 title "Z" with lines
set terminal pngcairo
set output "graphYZ.png"
set grid
set title "Evolucao de Y e Z em t"
set xlabel "t"
set ylabel "Temperatura" 
replot
set terminal wxt
set output
splot [0:1][0:1][0:1] "output_rk4.dat" using 3:4:5 title "Trajetoria" with lines
set terminal pngcairo
set output "graphXXZ.png"
set title "Trajetoria do sistema no espaco de fases"
set grid
set xlabel "X"
set ylabel "Y" 
set zlabel "Z" 
replot
set terminal wxt
set output
