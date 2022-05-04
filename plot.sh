#!/bin/bash

# the columns of the file that are to be plotted
X_COLUMN=2
Y_COLUMN=3

# check that the number of arguments is correct
if [ $# -ne 2 ];
then
    echo "The following 2 arguments are required, in this order:"
    echo "1. the parameter the file contains(eg. \"position\" or \"acceleration\")"
    echo "2. the value of theta used(eg. -1.57)"
    exit 2
fi

# get the arguments
PARAMETER=$1
THETA=$2

# factor to scale by
SCALE_FACTOR=50000000000

# gnuplot stuff
gnuplot -persist <<-EOFMarker
	set title "Path of sail in the XY plane"
	set key bottom right
	set xlabel "x position"
	set ylabel "y position"
	set xzeroaxis linetype -1
	set yzeroaxis linetype -1
	set grid back
	plot "${PARAMETER}_NFP_${THETA}.txt" using ${X_COLUMN}:${Y_COLUMN}, "jupiter_position_${THETA}.txt" using ${X_COLUMN}:${Y_COLUMN} with lines
	set xrange [GPVAL_DATA_X_MIN-$SCALE_FACTOR:GPVAL_DATA_X_MAX+$SCALE_FACTOR]
	set yrange [GPVAL_DATA_Y_MIN-$SCALE_FACTOR:GPVAL_DATA_Y_MAX+$SCALE_FACTOR]
	replot
EOFMarker
