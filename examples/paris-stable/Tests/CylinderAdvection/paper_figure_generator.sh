#! /bin/bash
#set -x

LIMITE=38
START=0
DELTA=100

/bin/rm -f *.png

for ((a=START; a <= LIMITE ; a++)) # Doubles parenthèses, et "LIMITE" sans "$".
do

let frame=$DELTA*$a
if [ $a == 0 ]; then
    frame=000
fi
if [ $a -lt 10 ]; then
    frame='00'$frame
else
    frame='0'$frame
fi
cp "out/CVoF-00000-$frame.txt" toplot$a.tmp
done
gnuplot < paperframe.gp
cp tmp.png cylindertest.png

