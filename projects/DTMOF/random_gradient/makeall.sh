cd extreme
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 1' | cat - Makefile > temp && mv temp Makefile
make
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 2' | cat - Makefile > temp && mv temp Makefile
make
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 3' | cat - Makefile > temp && mv temp Makefile
make
make clean
cd ..
cd uniform
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 1' | cat - Makefile > temp && mv temp Makefile
make
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 2' | cat - Makefile > temp && mv temp Makefile
make
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 3' | cat - Makefile > temp && mv temp Makefile
make
make clean
cd ..
cd exponential
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 1' | cat - Makefile > temp && mv temp Makefile
make
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 2' | cat - Makefile > temp && mv temp Makefile
make
make clean
sed -i '1d' Makefile
echo 'TEST_NUM = 3' | cat - Makefile > temp && mv temp Makefile
make
make clean
cd ..
