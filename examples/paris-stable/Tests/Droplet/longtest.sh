#!/bin/bash
#set -x

./run_one_test.sh 32 2 0.022 2d-5 F

gnuplot < plot.gp
mv droplet.png ../Testreport

