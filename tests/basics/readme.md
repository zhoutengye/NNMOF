make single file with 
	make

make all tests with
	bash makeall.sh

Note: the every test may use different HDF5 input file. See the ipynb file

test sctipts:
	mpirun -np 4 --oversubscribe ./TEST-1 
	mpirun -np 4 --oversubscribe ./TEST-2 test2
	mpirun -np 4 --oversubscribe ./TEST-3 test3
	mpirun -np 4 --oversubscribe ./TEST-4 test4

Some data have to be generated with the **mpi_test.ipynb**. 

For test4, it uses same data as test3, however, one may have to rename **test3.h5** to **test4.h5**.
