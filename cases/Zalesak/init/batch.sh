source ~/HPC/.intel.sh  
ulimit -s unlimited
cp ../../inputs/input50.namelist .
cp ../../inputs/input100.namelist .
./Zalesak-init.exec input50
cp output.h5 ../../inputs/input50.h5
./Zalesak-init.exec input100
cp output.h5 ../../inputs/input100.h5