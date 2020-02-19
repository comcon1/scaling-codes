#!/usr/bin/env python3

from multiprocessing.pool import ThreadPool
import os
import time
import subprocess
import sys
import numpy as np
import time

of = False

def launch(filename, num, params):
    start_time = time.time()
    print('Started')
    lst = ['../%s' % filename, "--output-file", "simu_%d" % num] +\
      list(map(str,params))
    print(' '.join(lst))
    out = subprocess.Popen(lst, universal_newlines=True, stdout=subprocess.PIPE)
    z = out.wait()
    print('Return: %d' % z)
    print('Ended: %f' % (time.time() - start_time) )

# main body starter

if __name__ == '__main__':
    usage_string = 'Usage: ' + sys.argv[0] + ' program_file parameters_file [num_of_cores]'
    if len(sys.argv) < 4 or len(sys.argv) > 5:
        raise ValueError('Invalid number of command line arguments. ' + usage_string)


    program_file = sys.argv[1]
    params_file = sys.argv[2]
    try:
        cores = int(sys.argv[3])
    except (IndexError) as e:
        cores = 1
    
    fd = open(params_file, 'r')
    params_list = fd.readlines()
    fd.close()
    
    print('Paramer sets: %d', len(params_list))
    directory = 'Sim_' + params_file  + '_' + time.strftime('%d.%m.%y-%H.%M.%S')

    if not os.path.exists(directory):
        os.makedirs(directory)
    os.chdir(directory)

    tp = ThreadPool(cores)
    for i in range(len(params_list)):
        tp.apply_async(launch, [program_file, i, params_list[i].split() ])

    tp.close()
    tp.join()
