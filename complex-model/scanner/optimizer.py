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
#    smt   cfa  cfx  cft
# =======================

initVector = [ 
    0.000990214267170159,9.318719739708396e-05,4.542275518749989,10.385969331018355,
    0.21539389927593255,6.062205518675469,0.00088947580412008,3.8368121167519265,
    2.9935202002665857,0.0008366572847518542,0.0003665298336643424,41.510057191938614,
    9.818065814677865,0.0004235577417342571,0.0009069921643147143,1.4152161322301913 ]

initBounds = ( (0,1e-3), (0,1e-3), (0,10), (6,20), 
               (0,5), (1,30), (5e-5,1e-3), (0.1,4),
               (0,10), (0,1e-3), (0,1e-3), (20,150), 
               (3.,10.), (0,1e-2), (0,1e-3), (0,10.) )
#               0      1        2       3      4       5        6       7     8      9
equiRads =   [ 725.,   551.,   698.,   538.,  670.,   649.,    637.,   791., 785., 770. ]
#               K   |  Hf   |  -N2  | Hf-N2 |  K   |  -M3  |  -M3-N2 |  K  | -C  | -C-N2
targetVals = np.array([
             [ 256.,   327.,   288.,   214.,   239.,  386.,    309.,   296., 358., 361. ],
             [ None,  None,    None,   None,  1374.,  989.,    None,   None, None, None ]
                  ])

flog = open('logoptim.log', 'w')

def cbck():
    global nice_kill
    print('''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        I WILL %s
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<''' % ('STOP' if nice_kill else 'CONTINUE') )
    return nice_kill

def _wrapper(x,y,z):
    global flog
    answer = launcher_1m.calculate_error(x,y,z, log=flog)
    if cbck():
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
except KeyboardInterrupt:
    print(">> Interrupted correctly! <<")
except:
    raise

flog.close()
