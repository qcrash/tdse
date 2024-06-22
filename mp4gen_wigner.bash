#!/bin/bash
gnuplot -e "set yrange [-2:2]; set xrange [-2:2]; set term png size 1600,1200; do for [i=1:1000] {; set output 'wigner_'.i.'.png'; splot 'wigner_'.i nonuniform matrix with points}"
ffmpeg -f image2 -r 10.0 -i wigner_%d.png -pix_fmt yuv420p -qscale 1 out.mp4
