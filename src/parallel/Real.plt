
set term png
set output "real.png"
set xlabel "File number"
set ylabel "Real-Time"
set ydata time
set timefmt "%M:%S"
set yrange [0: ] 
plot "job-file" using 2 with points
