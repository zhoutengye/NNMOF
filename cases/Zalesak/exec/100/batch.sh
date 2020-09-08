cd MOFBFGS
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
cd MOFLemoineGN
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
cd MOFNN
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
cd MOFNNStab
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
cd PLIC
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
cd THINC
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
cd MOFBFGSNUMERICAL
make clean
make 
cp ../../../inputs/input100.namelist . 
cp ../../../inputs/input100.h5 . 
./Zalesak* input100
cd ..
