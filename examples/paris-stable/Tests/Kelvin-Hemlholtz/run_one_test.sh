#!/bin/bash
#set -x


if [ $# -lt 5 ]; then
    echo "missing arguments for Kelvin-Helmholtz Instability test"
    echo usage $0 ismono rho_g_over_rho_l wavenumber ini_wave_amp  nx 
    exit
fi

ismono=$1
rhog=$2
nwavenumber=$3
ini_wave_amp=$4
nx=$5


if [ $ismono == T ]; then
    setmono=mono
    echo "mono"
else
    setmono=parallel
fi

let ny=2*$nx; nz=2
npx=2; let npy=2*$npx; npz=1

if [ $setmono == mono ]; then
    npy=1; npz=1; npx=1
fi

zlength=`awk -v "anx=$nx" 'BEGIN {print 2./anx}'`

/bin/rm -fr out input

let npstart=$npx*$npy*$npz

sed s/NXTEMP/$nx/g testinput.template | sed s/NPXTEMP/$npx/g | sed s/NYTEMP/$ny/g | sed s/NPYTEMP/$npy/g | sed s/ZLENGTHTEMP/$zlength/g > testinput
sed s/RHO2TEMP/$rhog/g testinput > testinput-$nx-$rhog 
ln -s testinput-$nx-$rhog input

sed  s/IWATEMP/$ini_wave_amp/g testinputvof.template | sed s/NWAVENUMBER/$nwavenumber/g > inputvof 

mpirun --oversubscribe -np $npstart paris > tmpout 2>&1
cat tmpout >> tmpoutall
cp stats stats-$rhog
gnuplot < plotd.gp

exit 0
