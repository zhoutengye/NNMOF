#!/bin/bash
#set -x

rm -fR input out
ln -s inputshort input
mpirun -np 3 paris > tmpout
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`
rm -f input

cat stats | awk '{ print $1 , $13+$14 }' | tail +2 > test_Ek.txt
sh plot.sh

precision=2e-1
pariscompare test_Ek.txt reference_Ek.txt $precision 1 0

GREEN="\\033[1;32m"
NORMAL="\\033[0m"



