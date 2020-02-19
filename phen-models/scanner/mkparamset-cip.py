#!/usr/bin/python

dm = [0.1,1]
de = [1,5]
am1 = [1,10,20,100]
am2 = [0.1,0.01,1e-3]
am3 = [1e-4,5e-4,1e-3]
ae = [1e-4,1e-5,1e-3]
ae2 = [1e-1,1e-2,1e-3]
flx = [1,2,10]
be = [0.01,1e-3,1e-4]
t0 = [0.1,0.5,1]
h = [2,4]
L = [50,100]

fd = open('paramset-cip.txt', 'w')

for _dm in dm:
  for _de in de:
    for _am1 in am1:
      for _am2 in am2:
        for _am3 in am3:
          for _ae in ae:
            for _ae2 in ae2:
              for _flx in flx:
                for _be in be:
                  for _t0 in t0:
                    for _h in h:
                      for _L in L:
                        fd.write('%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %3d %4d 1e5\n' % \
                          (_dm, _de, _am1, _am2, _am3, _ae, _ae2, _flx, _be, _t0, _h, _L) )

fd.close()

