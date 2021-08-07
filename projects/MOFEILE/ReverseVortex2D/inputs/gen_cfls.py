import os
base_cfl = 0.5
grid_list = [16,32,64,128,256]
cfl_list = [0.01,0.05,0.1,0.15,0.2,0.25,0.5,0.8,1.0]

for cfl in cfl_list:
    dir_name = 'cfl_'+str(cfl)
    os.system('mkdir -p ' + dir_name)
    for grid in grid_list:
        nml = 'input' + str(grid) + '.namelist'
        f1 = open(nml,'r')
        f2 = open(dir_name + '/' + nml,'w')
        for line in f1:
            if 'dt =' in line:
                dt = float(line.split('=')[1])
                dt_new = dt / base_cfl * cfl
                line2 = 'dt = ' + str(dt_new) + '\n'
            elif 'output_step =' in line:
                step = int(line.split('=')[1])
                step_new = int(step * base_cfl / cfl)
                line2 = 'output_step = ' + str(step_new) + '\n'
            else:
                line2 = line
            f2.write(line2)
        f1.close()
        f2.close()
