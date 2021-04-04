#!/usr/bin/env python3

from multiprocessing.pool import ThreadPool
import os, time, subprocess, sys
import numpy as np
from numpy import pi
import os.path

of = False

def calculate_error(_params, _erads, _eres, log=None):
    '''
    _params -- vector of results
       0: na in control sample (Noggin2 production)
       1: nx in all samples (Noggin2 degradation in absence of MMP3)
       2: mmp3 amount in control sample
       3: sa in all samples (SMAD cascade intensity)
       4: sb in all samples (BMP amount for SMAD switch)
       5: smad_threshold determines bifurcation point for SMAD amount (can be
         correlated with *sa*)
       6: nu determines rate of delta reaction to other agents
       7: h determines the strength of delta hysteresis
       8: tb determines the strength of chordin degradation inhibition by mmp3.
        Less tb - more inhibition.
       9: ca, production rate of Chordin in all samples
      10: cx, degratation rate of Chordin in all samples (in absence of MMP3)
      11: s0 for definition of de_init parameter
      12: (psa-cth) delta of threshold determines value of SMAD activation, which defines the edge of
        the somitic mesoderm (CACT expr. area)
      13: efa, ENAF production rate
      14: efx, ENAF degradation rate
      15: eft, ENAF threshold for cell fate switching
    _erads -- vector of equisurf. radii
    _eres -- vector of results
    @RETURN: float
    '''

    T = 10

    pna = _params[0];      pnx = _params[1];  pm3 = _params[2];  psa = _params[3]
    psb = _params[4];      pst = _params[5];  pnu = _params[6];  ph  = _params[7]
    ptb = _params[8];      pca = _params[9];  pcx = _params[10]; s0 = _params[11]
    pct = psa-_params[12]; efa = _params[13]; efx = _params[14]; eft = _params[15] 

    '''
             T     L      na      nx   mmp3      sa   sb  smthr  nu    h   tb    ca     cx   dinit
             0   1        2       3     4        5    6     7     8    9   10    11     12    13 
wh.      0   T   pi*r    pna     pnx   pm3      psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
hf.      1   T   pi*r    pna     pnx   pm3/10.  psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
wh. -N   2   T   pi*r    pna/5.  pnx   pm3      psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
hf. -N   3   T   pi*r    pna/5.  pnx   pm3/10.  psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
wh.      4   T   pi*r    pna     pnx   pm3      psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
wh. -M   5   T   pi*r    pna     pnx   pm3/5.   psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
wh. -M-N 6   T   pi*r    pna/5.  pnx   pm3/5.   psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
wh.      7   T   pi*r    pna     pnx   pm3      psa  psb   pst   pnu   ph  ptb   pca    pcx   5*s0/r^2/pi
wh. -C   8   T   pi*r    pna     pnx   pm3      psa  psb   pst   pnu   ph  ptb   pca/5. pcx   5*s0/r^2/pi
wh. -C-N 9   T   pi*r    pna/5.  pnx   pm3      psa  psb   pst   pnu   ph  ptb   pca/5. pcx   5*s0/r^2/pi

 PCT is not forwarded to C++ code because it acts at the stage of result interpretation
    '''
    _r = _erads
    param_array = np.array([ \
        [   T, pi*_r[0] ,pna   , pnx , pm3    , psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[0]**2, efa, efx, eft ],\
        [   T, pi*_r[1] ,pna   , pnx , pm3/10., psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[1]**2, efa, efx, eft ],\
        [   T, pi*_r[2] ,pna/5., pnx , pm3    , psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[2]**2, efa, efx, eft ],\
        [   T, pi*_r[3] ,pna/5., pnx , pm3/10., psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[3]**2, efa, efx, eft ],\
        [   T, pi*_r[4] ,pna   , pnx , pm3    , psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[4]**2, efa, efx, eft ],\
        [   T, pi*_r[5] ,pna   , pnx , pm3/5. , psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[5]**2, efa, efx, eft ],\
        [   T, pi*_r[6] ,pna/5., pnx , pm3/5. , psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[6]**2, efa, efx, eft ],\
        [   T, pi*_r[7] ,pna   , pnx , pm3    , psa, psb,  pst,  pnu,  ph, ptb,  pca   ,  pcx, 1.59e-5*s0*_r[7]**2, efa, efx, eft ],\
        [   T, pi*_r[8] ,pna   , pnx , pm3    , psa, psb,  pst,  pnu,  ph, ptb,  pca/5.,  pcx, 1.59e-5*s0*_r[8]**2, efa, efx, eft ],\
        [   T, pi*_r[9] ,pna/5., pnx , pm3    , psa, psb,  pst,  pnu,  ph, ptb,  pca/5.,  pcx, 1.59e-5*s0*_r[9]**2, efa, efx, eft ],\
            ])
    nexp = param_array.shape[0]
    assert(nexp == _eres.shape[1])
    log.write('Parameters vector: %s\n' % (' '.join(map(str, _params)))  )
    log.write('Radii vector: %s\n' % (' '.join(map(str, _erads)))  )
    # param_strings - parameters ready for getting to apply_async 
    param_strings = [[] for x in range(nexp)]
    for i in range(nexp):
        param_strings[i] += [ "-O" ]
        param_strings[i] += [ "_current_/exp%d.txt" % (i) ]
        param_strings[i] += [ "--input-params" ]
        param_strings[i] += map(str, param_array[i,:])
        param_strings[i] += [ "--change-d" ]
        param_strings[i] += [ "7.5", "0.5" ]
    # check _current_ folder
    #print(param_strings)
    try:
        _fd = open('_current_/1', 'w')
        _fd.write('check')
        _fd.close()
    except IOError as e:
        print("Cant write files into _current_")
        print(e)
        sys.exit(1)

    try:
        pass
        perform_simulations(param_strings, log)
    except Exception as e:
        print("Error in performing simulations. Killing..")
        print(e)
        sys.exit(1)

    try:
        g_results = extract_results_delta(param_array[:,1], "_current_", log)
        s_results = extract_results_smad(param_array[:,1], "_current_", log, \
                psa, psb, pct)
    except Exception as e:
        print("Error in extracting results. Killing..")
        raise e

    # sum CHD ignoring None
    rSQ = 0; _n = 0
    for ires in range(nexp):
        if not (_eres[0,ires] is None):
            rSQ += (_eres[0,ires] - g_results[ires])**2
            _n += 1
    
    rVal = 3*rSQ/_n

    # sum CACT ignoring None
    rSQ = 0; _n = 0;
    for ires in range(nexp):
        if not (_eres[1,ires] is None):
            rSQ += (_eres[1,ires] - s_results[ires])**2
            _n += 1
    
    rVal += rSQ/_n
    print(' RESULT DELTA VEC: %s\n' % (' '.join(map(str, g_results)) ))
    print(' RESULT SMAD VEC: %s\n' % (' '.join(map(str, s_results)) ))
    log.write('Resulting Delta: %s\n' % \
            (' '.join(map(lambda x: '%5.1f' % x, g_results)))  )
    log.write('Resulting  Smad: %s\n' % \
            (' '.join(map(lambda x: '%5.1f' % x, s_results)))  )
    log.write('Summary result: %f\n' % (rVal))
    log.flush()

    return rVal

# performing single simulation series
# parralelizing it into ncores
# paramSL - list of strings
# log - opened descriptor for log writing
def perform_simulations(paramSL, log):
    program_file = '../solver/tw.x'  
    try:
        num_cores = int(os.environ['L1M_NTHREADS'])
    except:
        num_cores = 2
    print('Use %d number of cores' % (num_cores) )
    nexp = len(paramSL)

    tp = ThreadPool(num_cores)
    for i in range(nexp):
        # launch(program_file, params, repeats)
        # launch(program_file, i, params_list[i,:] )
        tp.apply_async(launch, ( program_file, paramSL[i] ) )
    tp.close()
    tp.join()
    pass


# Extracting the main result from the trajectories
# Ls - array of simulation lengths, L (mcm)
# dirname - folder with trajectory files
# log - log descriptor
def extract_results_delta(Ls, dirname, log):
    nexp = len(Ls)
    res = np.zeros(nexp)
    try:
        for i in range(nexp):
            try:
                a = np.loadtxt(os.path.join(dirname, 'exp%d.txt' % (i)))
            except ValueError as e:
                print('Can\'t load the matrix')
                print('THE SOLUTION DIVERGES!')
                # Do not break calculation in this case;
                # just put -666 for identification of the divergence
                res[i] = -666
                continue
            Nx = (a.shape[1]-1) / 5.
            assert(abs(Nx - int(Nx)) < 1e-8) 
            Nx = int(Nx)
            print('Matrix loaded: ', Nx)
            xx = np.linspace(0, Ls[i], Nx)
            de = a[:,1+3*Nx:1+4*Nx] # delta distribution
            grd = np.gradient(de, axis=1)
            if (grd[-1,:].max() < 0.01):
                print('Maximum gradient: ', grd[-1,:].max())
                res[i] = -666
                continue
            delta_borders = xx[np.argmax(grd, axis=1)]
            res[i] = delta_borders[-1]
    except IOError as e:
        print('File does not exist!')
        raise e

    return res

def extract_results_smad(Ls, dirname, log, s_a, s_b, ctsh):
    nexp = len(Ls)
    res = np.zeros(nexp)
    try:
        for i in range(nexp):
            try:
                a = np.loadtxt(os.path.join(dirname, 'exp%d.txt' % (i)))
            except ValueError as e:
                print('Can\'t load the matrix')
                print('THE SOLUTION DIVERGES!')
                # Do not break calculation in this case;
                # just put -666 for identification of the divergence
                res[i] = -666
                continue
            Nx = (a.shape[1]-1) / 5.
            assert(abs(Nx - int(Nx)) < 1e-8) 
            Nx = int(Nx)
            print('Matrix loaded: ', Nx)
            xx = np.linspace(0, Ls[i], Nx)
            bm = a[:,1+4*Nx:] # BMP distribution
            sm = s_a*(bm/s_b)**2/(1+ (bm/s_b)**2) # SMAD distribution
            cact_border_i = sm[-1,:].searchsorted(ctsh)
            res[i] = Ls[i] if cact_border_i == Nx else xx[cact_border_i]
    except IOError as e:
        print('File does not exist!')
        raise e

    return res


# Launching the process
# * progname - the full path to exec file
# * params - list of ALL the parameters
def launch(progname, params):
    start_time = time.time()
    print('Started..')
    lst = [ progname ] +\
      list(map(str,params))
    print(' '.join(lst))
    out = subprocess.Popen(lst, universal_newlines=True, stdout=subprocess.PIPE)
    z = out.wait()
    print('Return: %d' % z)
    # while True:
    #   row = out.stdout.readline()
    #   if row == '':
    #     break
    #   print(row,end='')
    print('Ended: %f' % (time.time() - start_time) )
    return


# main body starter
if __name__ == '__main__':

    tsV = np.array([
        [ 256., 327., 288., 214., 239., 386., 309., 296., 358., 361. ],
        [ None, None, None, None,1374., 989., None, None, None, None ]
                    ])

    rsE = np.array(
        [ 725., 551., 698., 538., 670., 649., 637., 791., 785., 770. ])

    # =======================
    #    pna  pnx  pm3  psa 
    #    psb  pst  pnu  ph  
    #    ptb  pca  pcx  pde
    #    pct  efa  efx  eft  
    # =======================
    initVector = np.array([
 0.0009825491410708322, 7.66590953454e-05, 4.224027984344428, 10.3581808483661,
 0.214155408367835, 5.6896492525618285, 0.0004998229601811062, 3.98329388118285,
 3.2188620413130344, 0.0008449489223685573, 0.00033669412739256114, 36.82257415654826,
 0.5313070938965305, 1.97336462453e-05, 6.75165317732e-05, 1.4323111044574848 ])
    
    er = calculate_error(initVector, rsE.tolist(), tsV, log=sys.stdout)
    print('Calculated error: ', er)


