wget https://hdf-wordpress-1.s3.amazonaws.com/wp-content/uploads/manual/HDF5/HDF5_1_12_0/source/hdf5-1.12.0.tar.gz
tar -xzvf hdf5-1.12.0.tar.gz
cd hdf5-1.12.0
./configure --enable-fortran --enable-parallel
make install
