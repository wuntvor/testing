#set noborder
#set nokey
#set style line 1 lw 4
set title "Performance Comparison of LBC and BEST"
set xlabel "Total nodes  log(lx*ly*lz)"
set ylabel "Performance MLUP/s"
plot "performance.res" using 1:3 with lines title "LBC", "../lbbench/utl/benchmark/results/performance/trats_performance.dat" using 1:12 with lines title "BEST"
set size 1.0, 0.6
set terminal postscript portrait enhanced color dashed lw 1 "Helvetica" 14
set output "my-plot.ps"
replot
set terminal x11
#replot
set size 1,1
#pause -1
