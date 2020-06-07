#! /bin/bash
#set -x
export LANG=en_EN


if [ $# -lt 2 ]; then
    echo "missing arguments"
    echo usage $0  nx rhol
    exit
fi

nx=$1
rhol=$2

zlength=`awk -v "anx=$nx" 'BEGIN {print 2./anx}'`


/bin/rm -fr out input
mv stats stats.old   
sed s/NXTEMP/$nx/g testinput.template | sed s/ZLENGTHTEMP/$zlength/g | sed s/RHO2TEMP/$rhol/g  > testinput-$nx-$rhol
#sed s/ZLENGTHTEMP/$zlength/g  testinput.template | sed s/RHO2TEMP/$rhog/g  > testinput
ln -s testinput-$nx-$rhol input

 
mpirun -np 4 paris > tmpout
