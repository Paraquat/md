#!/usr/local/bin/python3

import argparse
import re
import sys
import scipy as sp
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import frame
from pdb import set_trace

parser = argparse.ArgumentParser(description='Analyse a MD trajectory',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--dump', action='store', dest='dumpfile',
                    default='dump.out', help='MD dump file')
parser.add_argument('-s', '--skip', action='store', dest='nskip', default=0,
                    help='Number of equilibration steps to skip')

opts = parser.parse_args()

nframes = 0
with open(opts.dumpfile, 'r') as infile:
  while True:
    if not line:
      break
    line = infile.readline()
