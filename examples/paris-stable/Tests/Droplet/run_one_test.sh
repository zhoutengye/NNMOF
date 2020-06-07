#!/bin/bash
#set -x


if [ $# -lt 4 ]; then
    echo "missing arguments"
    echo usage $0 nx npx endtime dt [MomCons: T/F]
    exit
fi

nx=$1
npx=1
endtime=$3
dt=$4

ny=$nx
nz=$nx

npy=2
npz=2
nfilter=1
plotlabel=MomCons
momconstrue=T


if [ $# -gt 4 ]; then
    if [ $5 == F ]; then
	echo "using nonMomCons"
	momconstrue=F
	plotlabel=nonMomCons
    else 
	if [ $5 == T ]; then
	    echo "using MomCons"
	else
	    echo "$0: usage error for parameter 5 = $4"
	    exit 1
	fi
    fi
fi


sed s/MOMCONSTEMP/$momconstrue/g inputvof.template | sed s/NFILTERTEMP/$nfilter/g > inputvof 
sed s/MOMCONSTEMP/$plotlabel/g plot.gp.template > plot.gp


if [ "$HAVE_SILO" == 1 ]; then
  echo "we have silo"
  out=5
else
  echo "we do not have silo"
  out=$2
fi


rm -fR out input
sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NZTEMP/$nz/g  | sed s/NPZTEMP/$npz/g  | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g  | sed s/ENDTEMP/$endtime/g | sed s/DTTEMP/$dt/g > testinput.tmp
ln -s testinput.tmp input
let npstart=$npx*$npy*$npz


if [ $npstart -gt 4 ]; then
    mpirun --oversubscribe -np $npstart paris > tmpout # 2>&1
else
    mpirun -np $npstart paris > tmpout # 2>&1
fi

# mpirun --oversubscribe -np $npstart paris > tmpout # 2>&1
echo `awk ' /Step:/ { cpu = $8 } END { print "cpu = " cpu } ' < tmpout`

precision=2e-1
pariscompare out/droplet-test-vel.txt reference.txt $precision 1 0

