cd MOFBFGS
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd ELVIRA
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd MOFLemoineGN
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd MOFNN
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd MOFNNStab
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd PLIC
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd THINC
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
cd MOFBFGSNUMERICAL
rm Linear* 
make clean
make 
cp ../../inputs/input.namelist . 
cp ../../inputs/input.h5 . 
./Linear* input
cd ..
