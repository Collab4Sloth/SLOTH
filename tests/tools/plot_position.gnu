set datafile separator ","

set xlabel "Time [s]"
set ylabel "Position [m]"
set grid
set terminal pngcairo
set output "output.png"

plot "data.csv" using 1:4 skip 1 with lines
pause -1
