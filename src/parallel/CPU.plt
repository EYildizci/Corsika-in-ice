set term png
set output "cpu.png"
set xlabel "File number"
set ylabel "CPU-Time"

plot "job-file" using 4 with points
