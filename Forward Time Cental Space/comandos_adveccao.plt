set term png
set output 'AdvectionGraph.png'
set grid
set style data lines
unset key
set xlabel 'X'
set ylabel 'Time'
splot 'comandos_adveccao.plt' using 1:2:3 with lines
quit
