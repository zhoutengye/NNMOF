#! /bin/bash
#set -x

run=1
if [ $# -lt 4 ]; then
    echo "missing arguments"
    echo usage $0 ismono nwavenumber ini_wave_amp  nx \[run=1\]
    exit 1
fi

ismono=$1
nwavenumber=$2
ini_wave_amp=$3
nx=$4
if [ $# -eq 5 ]; then
run=$5
fi
momcons=`grep DoMOMCONS testinputvof.template | tr -d " "`
scheme=`cat testinputvof.template | awk -F "=" '/VOF_ADVECT/ {print $2}' | awk '{print $1}'`
scheme=`echo $scheme | tr -d "\'" | sed s/_/-/g`
if [ "$scheme" = "Dick-Yue" ] ; then
    scheme=Weymouth-Yue
fi
# ---- what is this for ? ----
ini_wave_amp_sp=`echo $ini_wave_amp | sed s/d/e/`

# amplitude constant for plotting exponential growth
# amplitude=$ini_wave_amp_sp
# ---- end ? --- 
amplitude=$ini_wave_amp

rhoglist=" 1    0.5   0.2   0.1    0.05"
cp multiplot.template.gp temp.gp
for rhog in $rhoglist; do 
if [ $run -ne 0 ]; then
    ./run_one_test.sh $ismono $rhog $nwavenumber $amplitude  $nx 
fi
    sed s/SRHOG/stats-$rhog/ temp.gp |  sed s/RHOG/$rhog/ > temp1.gp
    mv temp1.gp temp.gp
done

drhog=1
# Compute growth rate assuming ∆U = 1
stemp=`awk -v "r=$drhog" -v "w=$nwavenumber" 'BEGIN {ak=2.*w*atan2(0, -1); print sqrt(r)*ak}'`
sed s/ATEMP/$ini_wave_amp/ temp.gp |  sed s/STEMP/$stemp/ |  sed s/DRHOG/$drhog/ > temp1.gp
drhog=0.05
# Compute again growth rate 
stemp=`awk -v "r=$drhog" -v "w=$nwavenumber" 'BEGIN {ak=2.*w*atan2(0, -1); print sqrt(r)*ak}'`
sed s/ATEMP/$ini_wave_amp/ temp1.gp |  sed s/STEMP/$stemp/ |  sed s/DRHOG/$drhog/ | sed s/TITLE/KH' 'vmax' 'growth' 'for' 'nx=$nx' 'wvnr=$nwavenumber' '$momcons' '$scheme/ > temp.gp
mv temp.gp multiplot.gp
gnuplot < multiplot.gp

# -------- second set : low air density ----------

cp multiplot.template.gp temp.gp
rhoglist2="0.02 0.01  0.005 0.002  0.001"
for rhog in $rhoglist2; do 
if [ $run -ne 0 ] ; then 
    ./run_one_test.sh $ismono $rhog $nwavenumber $amplitude $nx 
fi 
    sed s/SRHOG/stats-$rhog/ temp.gp |  sed s/RHOG/$rhog/ > temp1.gp
    mv temp1.gp temp.gp
done
drhog=0.02
# Compute growth rate assuming ∆U = 1
stemp=`awk -v "r=$drhog" -v "w=$nwavenumber" 'BEGIN {ak=2.*w*atan2(0, -1); print sqrt(r)*ak}'`
sed s/ATEMP/$ini_wave_amp/ temp.gp |  sed s/STEMP/$stemp/ |  sed s/DRHOG/$drhog/ > temp1.gp
drhog=0.001
# Compute again growth rate 
stemp=`awk -v "r=$drhog" -v "w=$nwavenumber" 'BEGIN {ak=2.*w*atan2(0, -1); print sqrt(r)*ak}'`
sed s/ATEMP/$ini_wave_amp/ temp1.gp |  sed s/STEMP/$stemp/ |  sed s/DRHOG/$drhog/ | sed s/TITLE/KH' 'vmax' 'growth' 'for' 'nx=$nx' 'wvnr=$nwavenumber' '$momcons' '$scheme/ > temp.gp
mv temp.gp multiplot.gp
gnuplot < multiplot.gp
