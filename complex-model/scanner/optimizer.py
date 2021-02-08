#!/usr/bin/python3

from scipy.optimize import minimize
from scipy.optimize import dual_annealing
import launcher_1m
import numpy as np
import sys, signal

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

# =======================
#    pna   pnx  pm3  psa 
#    psb   pst  pnu  ph  
#    ptb   pca  pcx  pde 
#    smt   efa  efx  eft
# =======================

initVector = [
 0.0009825491410708322, 7.66590953454e-05, 4.224027984344428, 10.3581808483661,
 0.214155408367835, 5.6896492525618285, 0.0004998229601811062, 3.98329388118285,
 3.2188620413130344, 0.0008449489223685573, 0.00033669412739256114, 36.82257415654826,
 0.5313070938965305, 1.97336462453e-05, 6.75165317732e-05, 1.4323111044574848
]

initBounds = ( (0,1e-3), (0,1e-3), (0,10), (6,20), 
               (0,5), (1,30), (5e-5,1e-3), (0.1,4),
               (0,10), (0,1e-3), (0,1e-3), (20,150), 
               (0,5.), (0,1e-2), (0,1e-3), (0,10.) )
#               0      1        2       3      4       5        6       7     8      9
equiRads =   [ 725.,   551.,   698.,   538.,  670.,   649.,    637.,   791., 785., 770. ]
#               K   |  Hf   |  -N2  | Hf-N2 |  K   |  -M3  |  -M3-N2 |  K  | -C  | -C-N2
targetVals = np.array([
             [ 256.,   327.,   288.,   214.,   239.,  386.,    309.,   296., 358., 361. ],
             [ None,  None,    None,   None,  1374.,  989.,    None,   None, None, None ]
                  ])

flog = open('logoptim_mc2.log', 'w')

def cbck(param_vector, f_res, ctx):
    global nice_kill
    global flog
    print('''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       STEP FINISHED I WILL %s
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<''' % ('STOP' if nice_kill else 'CONTINUE') )
    flog.write("\n>> STEP RESULTS <<\n\n")
    flog.write("Found params: \n %s,\n %s,\n %s,\n %s\n" %
      ( ', '.join(map(str,param_vector[0:4])) ,  
        ', '.join(map(str,param_vector[4:8])),  
        ', '.join(map(str,param_vector[8:12])), 
        ', '.join(map(str,param_vector[12:16])) ) )

    flog.write("Found error: %-10.3f" % (f_res) )

    flog.write("\n>> STEP RESULTS <<\n")

    flog.flush()
    return nice_kill

def _wrapper(x,y,z):
    global flog
    global nice_kill
    # for imitational run
    # answer = np.random.random()*0.01 + np.sum(np.array(x)**2)
    answer = launcher_1m.calculate_error(x,y,z, log=flog)
    if nice_kill:
        raise KeyboardInterrupt('INTERRUPT')
    return answer

try:
    res = dual_annealing(_wrapper, \
        x0=initVector, \
        args=( equiRads, targetVals, ), \
        bounds=initBounds,
        initial_temp=100, \
        visit=1.5, \
        callback=cbck \
        )
    nice_kill = True
    flog.write("\n === FINAL STEP INFORMATION === \n")
    cbck(res.x, res.fun, 0)
    res_str = res.message
    print('\n --> '.join( res_str) )
    flog.write('\n'.join(res_str))
except KeyboardInterrupt:
    print(">> Interrupted correctly! <<")
except:
    raise
finally:
    flog.close()
