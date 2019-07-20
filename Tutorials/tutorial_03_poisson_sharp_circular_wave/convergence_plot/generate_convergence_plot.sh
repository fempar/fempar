#!/bin/bash
# !!!NOTE: This bash shell script requires gnuplot and epstopdf

echo "set term postscript eps enhanced color" > plot.plot
echo "set output \"uniform_mesh_error.eps\"" >> plot.plot
echo "set grid" >> plot.plot
echo "set yrange [1e-06:100]" >> plot.plot
echo "set xrange [0:1000]" >> plot.plot
echo "set key bottom left" >> plot.plot
echo "set logscale y" >> plot.plot
echo "set xtics 0,100,1000" >> plot.plot
echo "set ytics format \"%.1e\"" >> plot.plot
echo "set ylabel \"True error in energy norm (H1 seminorm)\"" >> plot.plot
echo "set xlabel \"Square root of the number of DOFs\"" >> plot.plot
echo "set title \"2D Poisson sharp circular wave front. CG, ALPHA=200, RADIUS=0.7, CENTER=(-0.05,-0.05)\"" >> plot.plot
plot_command="plot 'uniform_mesh_error_order_1' title 'order 1 uniform' with linespoints ps 1.8 lt 1 dt 1 lw 3 pt 9 lc 2, \
	           'uniform_mesh_error_order_2' title 'order 2 uniform' with linespoints ps 1.6 lt 1 dt 1 lw 3 pt 5 lc 3,  \
	           'uniform_mesh_error_order_4' title 'order 4 uniform' with linespoints ps 1.8 lt 1 dt 1 lw 3 pt 13 lc 4,  \
	           'uniform_mesh_error_order_8' title 'order 8 uniform' with linespoints ps 1.6 lt 1 dt 1 lw 3 pt 7 lc 5" 
		    
echo $plot_command >> plot.plot

gnuplot plot.plot
epstopdf uniform_mesh_error.eps
rm -f uniform_mesh_error.eps
rm -f plot.plot
