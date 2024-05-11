#!/bin/bash
gnuplot -e "set yrange [-2:2]; set zrange [-2:2]; set arrow from -1,0,0
to 1,0,0; set term png size 1600,1200; do for [i=1:100] {; set output
'psi_'.i.'.png'; splot 'psi_'.i using 1:2:3 w impulses w l}"
ffmpeg -f image2 -r 10.0 -i psi_%d.png -pix_fmt yuv420p -qscale 1 out.mp4
