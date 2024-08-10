#!/bin/bash
gnuplot -e "set yrange [-2:2]; set xrange [-2:2]; set zrange [0:0.3];
set xlabel 'x'; set ylabel 'p'; set zlabel 'Prob. Density'; set
xzeroaxis; set yzeroaxis; set contour; set cntrparam levels 11; set
cntrparam levels incremental 0,0.03,0.3; set term png size 1600,1200; do for [i=1:1000] {; set output 'sb_'.i.'.png'; splot 'sb_'.i nonuniform matrix with points}"
ffmpeg -f image2 -r 10.0 -i sb_%d.png -pix_fmt yuv420p -qscale 1 out.mp4
