plot "output_thomas_analytic.dat" using 2:3 title "Analitica" with lines
set terminal pngcairo
set output "Temperatura1D_analytic.png"
set grid
set title "Temperatura ao longo da placa"
set xlabel "X"
set ylabel "Temperatura"
replot
set terminal wxt
set output
plot "output_thomas_crank_nicolson.dat" using 2:3 title "Crank-Nicolson" with lines
set terminal pngcairo
set output "Temperatura1D_numeri.png"
set grid
set title "Temperatura ao longo da placa"
set xlabel "X"
set ylabel "Temperatura"
replot
set terminal wxt
set output