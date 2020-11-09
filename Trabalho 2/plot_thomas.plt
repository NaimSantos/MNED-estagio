plot "output_thomas_analytic.dat" using 2:3 title "Solucao Analitica" with lines
set terminal pngcairo
set output "Temperatura1D.png"
set grid
set title "Temperatura ao longo da placa"
set xlabel "X"
set ylabel "Temperatura" 
replot
set terminal wxt
set output