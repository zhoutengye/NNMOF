cd MOFBFGS
make clean
make 
cp ../../../inputs/input50.namelist . 
cp ../../../inputs/input50.h5 . 
./Zalesak* input50
cd ..
cd MOFLemoineGN
make clean
make 
cp ../../../inputs/input50.namelist . 
cp ../../../inputs/input50.h5 . 
./Zalesak* input50
cd ..
