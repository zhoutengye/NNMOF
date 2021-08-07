case_list = [16,32,64,128,256]

f = open('run_all.sh','w')
for ca in case_list:
    case = str(ca)
    f.write('cp ../inputs/input' + case + '.namelist .\n')
f.write('\n')
for ca in case_list:
    case = str(ca)
    f.write('./RV-init.exec input' + case + '\n')
    f.write('cp output.h5 ../inputs/input' + case + '.h5\n')

f.close()
