import os

BENCH_case = 'Benchmark'
cfl_list = ['0.5','0.8','1.0']
grid_list = ['16','32','64','128']

for cfl in cfl_list:
    cfl_dir = 'cfl_' + str(cfl)
    os.system('mkdir -p ' + cfl_dir)
    for case in grid_list:
        case_dir = cfl_dir + '/' + case
        print('cp ' + BENCH_case + '/gen_cases.py ' + case_dir)
        os.system('mkdir -p ' + case_dir)
        os.system('cp ' + BENCH_case + '/gen_cases.py ' + case_dir + '/')
        os.system('cp -r ' + BENCH_case + '/setup.json ' + case_dir + '/')
        f = open(case_dir+'/gen_cases.py','r')
        new_lines = []
        for line in f:
            if 'reso =' in line:
                new_lines.append("reso = '" + case + "'\n")
            elif 'cfl =' in line:
                new_lines.append("cfl = '" + cfl + "'\n")
            else:
                new_lines.append(line)
        f.close()
        f2 = open(case_dir+'/gen_cases.py','w')
        for line in new_lines:
            f2.write(line)
        f2.close()

f0 = open('run_cfl050810.sh','w')
f01 = open('clean_cfl050810.sh','w')
for cfl in cfl_list:
    cfl_dir = 'cfl_' + str(cfl)
    f0.write('cd ' + cfl_dir + '\n')
    f0.write('bash run_all.sh' + '\n')
    f0.write('cd .. \n')
    f01.write('rm -r ' + cfl_dir + '\n')
    f = open(cfl_dir + '/run_all.sh','w')
    f2 = open(cfl_dir + '/clean_all.sh','w')
    for case in grid_list:
        f.write('cd ' + case + '\n')
        f.write('python gen_cases.py\n')
        f.write('bash run_cases.sh\n')
        f.write('cd ..\n')
        f2.write('rm -r ' + case + '\n')
    f.close()
    f2.close()
