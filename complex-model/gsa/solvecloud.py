#!/usr/bin/python3

import numpy as np
import sys, signal
sys.path.append('../scanner')
import launcher_1m

nice_kill = False

def exit_gracefully(signum, frame):
    global nice_kill
    print('''
============================================
PROCESS WILL BE TERMINATED AT THE NEXT STEP
============================================
''')
    nice_kill = True

signal.signal(signal.SIGINT, exit_gracefully)
signal.signal(signal.SIGTERM, exit_gracefully)

#               0      1        2       3      4       5        6       7     8      9
equiRads =   [ 725.,   551.,   698.,   538.,  670.,   649.,    637.,   791., 785., 770. ]
#               K   |  Hf   |  -N2  | Hf-N2 |  K   |  -M3  |  -M3-N2 |  K  | -C  | -C-N2
targetVals = np.array([
             [ 256.,   327.,   288.,   214.,   239.,  386.,    309.,   296., 358., 361. ],
             [ None,  None,    None,   None,  1374.,  989.,    None,   None, None, None ]
                  ])

flog = open('saltelli-5-100.log', 'w')
fval = open('values_sobol-5-100.txt', 'w')

N_START = 0
cfd  = open('saltelli-5-100.csv', 'r')

def _wrapper(x):
    global flog
    global nice_kill
    # for imitational run
    # answer = np.random.random()*0.001 + np.sum(np.array(x)**2); print('.', end='')
    answer = launcher_1m.calculate_error(x, equiRads, targetVals, log=flog)
    if nice_kill:
        raise KeyboardInterrupt('INTERRUPT')
    return answer

try:
    calcFlag = False
    while True:
        line = cfd.readline()
        if line == '':
            print('End-of-file was reached!')
            break
        splLine = line.split(',')
        if splLine[0] == str(N_START):
            print('Line found!')
            calcFlag = True
        if not calcFlag:
            continue
        curLine = int(splLine[0])
        params = list(map(float,splLine[1:17]))
        curVal = _wrapper(params)
        fval.write('%d,%15.8e\n' % (curLine, curVal) )
        fval.flush()
except KeyboardInterrupt:
    print(">> Interrupted correctly! <<")
except:
    raise
finally:
    flog.close()
    cfd.close()
    fval.close()

