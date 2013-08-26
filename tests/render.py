#! /usr/bin/env python
#

import sys

import pylab

if __name__ == '__main__':
  file = sys.argv[1]
  pos = []
  forces = []
  execfile(file)
  print pos
  print forces

  X = [p[0] for p in pos]
  Y = [p[1] for p in pos]

  pylab.plot(X, Y)
  for f in forces:
    FX = [f[0][0], f[1][0]/100.]
    FY = [f[0][2], f[1][2]/100.]
    pylab.plot(FX, FY)

  pylab.show()
