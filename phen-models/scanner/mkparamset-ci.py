#!/usr/bin/python

dm = [10,20,30]
de = [10,20,30]
am1 = [1e-5,1e-4,1e-3,1e-2,.1,1]
am2 = [1e-3,1e-2,1e-1,1,10,1e2]
ae = [1e-5,1e-4,1e-3,1e-2,.1,1]
flx = [.1,1,10]
be = [1e-4,1e-3,1e-2]
t0 = [1e-5,1e-4,1e-3]
h = [2,4]
L = [50,100]

fd = open('paramset-ci.txt', 'w')

for _dm in dm:
  for _de in de:
    for _am1 in am1:
      for _am2 in am2:
        for _ae in ae:
          for _flx in flx:
            for _be in be:
              for _t0 in t0:
                for _h in h:
                  for _L in L:
                    fd.write('%10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f %3d %4d 1e5\n' % \
                      (_dm, _de, _am1, _am2, _ae, _flx, _be, _t0, _h, _L) )

fd.close()
