cp ../inputs/input16.namelist .
cp ../inputs/input32.namelist .
cp ../inputs/input64.namelist .
cp ../inputs/input128.namelist .
cp ../inputs/input256.namelist .

./RV-init.exec input16
cp output.h5 ../inputs/input16.h5
./RV-init.exec input32
cp output.h5 ../inputs/input32.h5
./RV-init.exec input64
cp output.h5 ../inputs/input64.h5
./RV-init.exec input128
cp output.h5 ../inputs/input128.h5
./RV-init.exec input256
cp output.h5 ../inputs/input256.h5
