# Each segment, delete the first line of makefile and append a 
# new line at the begining of the file
#
sed -i '1d' Makefile
echo 'TEST_NUM = 1' | cat - Makefile > temp && mv temp Makefile
make
make clean

sed -i '1d' Makefile
echo 'TEST_NUM = 2' | cat - Makefile > temp && mv temp Makefile
make
make clean
