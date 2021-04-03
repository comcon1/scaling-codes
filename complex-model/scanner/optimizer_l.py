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

# ===========================
#   pna       pnx  pm3  psa 
#   psb       pst  pnu  ph  
#   ptb       pca  pcx  pde 
#   psa-smt   efa  efx  eft
# ===========================

initVector = [ 
 0.0009907912633801976, 7.314768997741207e-05, 4.312724643922131, 10.366895759907635, 
 0.2145836479258079, 5.657215083250048, 0.0005588656869399085, 3.983260474504072, 
 3.185733880385387, 0.0008566201401132153, 0.000342246930252877, 37.19903211401452, 
 0.5400000000000, 1.977310697154443e-05, 6.760879769684623e-05, 1.429242155130801 ]

epsVector = np.array(initVector)/20.

initBounds = ( (0,1e-3), (0,1e-3), (0,10), (6,20), 
               (0,5), (1,30), (5e-5,1e-3), (0.1,4),
               (0,10), (0,1e-3), (0,1e-3), (20,150), 
               (0,5), (0,1e-2), (0,1e-3), (0,10.) )
#               0      1        2       3      4       5        6       7     8      9
equiRads =   [ 725.,   551.,   698.,   538.,  670.,   649.,    637.,   791., 785., 770. ]
#               K   |  Hf   |  -N2  | Hf-N2 |  K   |  -M3  |  -M3-N2 |  K  | -C  | -C-N2
targetVals = np.array([
             [ 256.,   327.,   288.,   214.,   239.,  386.,    309.,   296., 358., 361. ],
             [ None,  None,    None,   None,  1374.,  989.,    None,   None, None, None ]
                  ])

flog = open('logoptim.log', 'w')

def cbck(param_vector):
    global nice_kill
    global flog
    print('''
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
       STEP FINISHED I WILL %s
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<''' % ('STOP' if nice_kill else 'CONTINUE') )
    flog.write("\n>> STEP RESULTS <<\n\n")
    flog.write("Found params: \n %s\n %s\n %s\n %s\n" %
      ( ', '.join(map(str,param_vector[0:4])) ,  
        ', '.join(map(str,param_vector[4:8])),  
        ', '.join(map(str,param_vector[8:12])), 
        ', '.join(map(str,param_vector[12:16])) ) )
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
    res = minimize(_wrapper,
                   x0=initVector,
                   args=( equiRads, targetVals, ),
                   bounds=initBounds,
                   method='L-BFGS-B', 
                   callback=cbck, 
                   options = { 'maxcor': 30, 
                               'maxfun': 1e6,
                               'eps': epsVector,
                               'maxiter': 1e5,
                               'gtol': 1e-7,
                               'ftol': 1e-12 }
                  )
    nice_kill = True
    cbck(res.x)
    res_str = res.message.decode('utf-8')
    flog.write(res_str)
    print(res_str)
except KeyboardInterrupt:
    print(">> Interrupted correctly! <<")
except:
    raise

flog.close()
