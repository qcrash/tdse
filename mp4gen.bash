#!/bin/bash
gnuplot -e "set yrange [-2:2]; set xrange [-2:2]; set term png size 1600,1200; do for [i=1:1000] {; set output 'psi_'.i.'.png'; p 'psi_'.i using 2:3 w l}"
ffmpeg -f image2 -r 10.0 -i psi_%d.png -pix_fmt yuv420p -qscale 1 out.mp4
