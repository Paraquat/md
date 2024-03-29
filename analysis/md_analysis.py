#!/usr/local/bin/python3

import argparse
import sys
import re
import scipy as sp
import matplotlib.pyplot as plt
from frame import Frame
from md_tools import Pairdist, MSER, VACF, MSD
from pdb import set_trace

ha2ev = 27.211399

# Regular expressions
frame_re = re.compile('frame')
endframe_re = re.compile('end frame')
specblock_re = re.compile(r'%block ChemicalSpeciesLabel\n(.*?)\n%endblock',
                          re.M | re.S | re.I)

input_file = 'md.in'
coord_file = 'cell.in'

# Parsing functions

def strip_comments(line, separator):
  for s in separator:
    i = line.find(s)
    if i >= 0:
      line = line[:i]
  return line.strip()


def parse_input(input_file):
  input_params = {}
  with open(input_file, 'r') as input:
    for line in input:
      stripped = strip_comments(line, "#%!")
      if stripped:
        bits = stripped.split()
        if len(bits[1:]) == 1:
          input_params[bits[0]] = bits[1:][0]
        else:
          input_params[bits[0]] = bits[1:]

  return input_params

def parse_cell(cell_file):
  data = {}
  cell_re = re.compile("%BLOCK lattice_cart(.*?)%ENDBLOCK lattice_cart", re.S | re.M | re.I)
  pos_re = re.compile("%BLOCK positions_frac(.*?)%ENDBLOCK positions_frac", re.S | re.M | re.I)
  with open(cell_file, 'r') as infile:
    text = infile.read()
  lattext = re.findall(cell_re, text)[0].strip().splitlines()
  postext = re.findall(pos_re, text)[0].strip().splitlines()
  latvec = []
  for line in lattext:
    lat = [float(bit) for bit in line.split()]
    latvec.append(lat)
  data['latvec'] = sp.array(latvec)
  species = []
  coords = []
  for line in postext:
    s, x, y, z = line.split()
    species.append(s)
    coords.append([float(x), float(y), float(z)])
  natoms = len(coords)
  data['species'] = sp.array(species)
  data['coords'] = sp.array(coords)
  scount = {}
  for i in range(natoms):
    if data['species'][i] in scount.keys():
      scount[data['species'][i]] += 1
    else:
      scount[data['species'][i]] = 1
  data['species_count'] = scount
  data['nspecies'] = len(scount.keys())

def read_stats(stats_file, nstop):
  nstep = 0
  data = {}
  header = True
  with open(stats_file, 'r') as statfile:
    for line in statfile:
      if nstop != -1:
        if nstep > nstop:
          break
      if header:
        col_id = line.strip().split()
        for col in col_id:
          data[col] = []
        header = False
      else:
        bits = line.strip().split()
        for i, bit in enumerate(bits):
          if i==0:
            info = int(bit)
          else:
            info = float(bit)
          data[col_id[i]].append(info)
      nstep += 1
  return nstep, data

# Command line arguments
parser = argparse.ArgumentParser(description='Analyse a MD trajectory',
                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--dirs', nargs='+', default='.', dest='dirs',
                    action='store', help='Directories to analyse')
parser.add_argument('-f', '--frames', action='store', dest='framesfile',
                    default='Frames', help='MD frames file')
parser.add_argument('-s', '--stats', action='store', dest='statfile',
                    default='Stats', help='MD statistics file')
parser.add_argument('--skip', action='store', dest='nskip', default=0,
                    type=int, help='Number of equilibration steps to skip')
parser.add_argument('--stop', action='store', dest='nstop', default=-1, 
                    type=int, help='Number of last frame in analysis')
parser.add_argument('--equil', action='store', dest='nequil', default=0, 
                    type=int, help='Number of equilibration steps')
parser.add_argument('--vacf', action='store_true', dest='vacf', 
                    help='Plot velocity autocorrelation function')
parser.add_argument('--msd', action='store_true', dest='msd', 
                    help='Plot mean squared deviation')
parser.add_argument('--rdf', action='store_true', dest='rdf', 
                    help='Plot radial distribution function')
parser.add_argument('--stress', action='store_true', dest='stress', 
                    help='Plot stress')
parser.add_argument('--landscape', action='store_true', dest='landscape', 
                    help='Generate plot with landscape orientation')
parser.add_argument('--pub', action='store_true', dest='pub', 
                    help='Publication text size')
parser.add_argument('--nbins', action='store', dest='nbins', default=100,
                    help='Number of histogram bins')
parser.add_argument('--rdfwidth', action='store', dest='rdfwidth',
                    default=0.05, help='RDF histogram bin width (A)')
parser.add_argument('--rdfcut', action='store', dest='rdfcut', default=10.0,
                    help='Distance cutoff for RDF')
parser.add_argument('--dump', action='store_true', dest='dump', 
                    help='Dump secondary data used to generate plots')
parser.add_argument('--mser', action='store', dest='mser_var', default=None,
                    type=str, help='Compute MSER for the given property')

opts = parser.parse_args()
if (opts.vacf or opts.msd or opts.stress or opts.rdf):
  read_frames = True
else:
  read_frames = False

if opts.nequil == 0:
  opts.nequil = opts.nskip

# Parse the md.in parameters file
input_params = parse_input(input_file)
init_config = parse_cell(coord_file)
natoms = int(input_params['natoms'])
dt = float(input_params['dt'])

# Parse the statistics file
nsteps, data = read_stats(opts.statfile,opts.nstop)
avg = {}
std = {}
for key in data:
  data[key] = sp.array(data[key])
  avg[key] = sp.mean(data[key][opts.nequil:-1])
  std[key] = sp.std(data[key][opts.nequil:-1])
time = [float(s)*dt for s in data['step']]
data['time'] = sp.array(time)

# Plot the statistics
if opts.landscape:
  fig1, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2, sharex=True, figsize=(11,7))
  plt.tight_layout(pad=6.5)
else:
  fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True, figsize=(7,10))

ax1.plot(data['time'][opts.nskip:], data['pe'][opts.nskip:], 'r-', label='Potential energy')
ax1a = ax1.twinx()
ax1a.plot(data['time'][opts.nskip:], data['ke'][opts.nskip:], 'b-', label='Kinetic energy')
if input_params['ensemble'][2] == 't':
  if input_params['thermo_type'] == 'nhc':
    ax1a.plot(data['time'][opts.nskip:], data['nhc'][opts.nskip:], 'g-', label='NHC energy')
if input_params['ensemble'][1] == 'p':
  if 'mttk' in input_params['baro_type']:
    ax1a.plot(data['time'][opts.nskip:], data['box'][opts.nskip:], 'c-', label='Box energy')
  ax1a.plot(data['time'][opts.nskip:], data['pV'][opts.nskip:], 'm-', label='pV')
ax2.plot(data['time'][opts.nskip:], data['H\''][opts.nskip:])
ax2.plot((opts.nskip,data['time'][-1]), (avg['H\''],avg['H\'']), '-',
      label=r'$\langle H\' \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['H\''], std['H\'']))
ax3.plot(data['time'][opts.nskip:], data['T'][opts.nskip:])
ax3.plot((opts.nskip,data['time'][-1]), (avg['T'],avg['T']), '-',
      label=r'$\langle T \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['T'], std['T']))
ax4.plot(data['time'][opts.nskip:], data['P'][opts.nskip:], 'b-')
ax4.plot((opts.nskip,data['time'][-1]), (avg['P'],avg['P']), 'b--',
      label=r'$\langle P \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['P'], std['P']))
if input_params['ensemble'][1] == 'p':
  ax4a = ax4.twinx()
  ax4a.plot(data['time'][opts.nskip:], data['V'][opts.nskip:], 'r-')
  ax4a.plot((opts.nskip,data['time'][-1]), (avg['V'],avg['V']), 'r--',
        label=r'$\langle V \rangle$ = {0:>12.4f} $\pm$ {1:<12.4f}'.format(avg['V'], std['V']))
ax1.set_ylabel("E (Ha)")
ax2.set_ylabel("H$'$ (Ha)")
ax3.set_ylabel("T (K)")
ax4.set_ylabel("P (GPa)", color='b')
if input_params['ensemble'][1] == 'p':
  ax4a.set_ylabel("V ($a_0^3$)", color='r')
ax4.set_xlabel("time (fs)")
ax1.legend(loc="upper left")
ax1a.legend(loc="lower right")
ax2.legend()
ax3.legend()
ax4.legend(loc="upper left")
if input_params['ensemble'][1] == 'p':
  ax4a.legend(loc="lower right")
plt.xlim((opts.nskip,data['time'][-1]))
fig1.subplots_adjust(hspace=0)
fig1.savefig("stats.pdf", bbox_inches='tight')

# Plot MSER
if opts.mser_var:
  traj = MSER(nsteps, opts.mser_var, data[opts.mser_var])
  traj.get_mser()
  traj.plot_mser(data['step'])
  if opts.dump:
    traj.dump_mser(data['step'])

# Parse the frames file
if read_frames:
  nframes = 0
  newframe = True
  buf = ""
  time = []
  stress = []
  lat = []
  first_frame = True
  with open(opts.framesfile, 'r') as framesfile:
    while True:
      line = framesfile.readline()
      if not line:
        break
      if re.match(frame_re, line):
        n = int(line.split()[1])

      if re.match(endframe_re, line):
        newframe = True
        if n <= opts.nskip:
          continue
        else:
          nframes += 1
        if opts.nstop != -1:
          if n > opts.nstop:
            break
        sys.stdout.write("Processing frame {}\r".format(n))
        if first_frame:
          first_frame = False
          f1 = Frame(natoms,n)
          f1.parse_frame(buf)
          if opts.rdf:
            pairdist = Pairdist(natoms, init_config['nspecies'], opts.rdfcut, 
                                opts.rdfwidth, input_params['species'],
                                init_config['species_count'])
          if opts.vacf:
            c = VACF(natoms, dt, f1)
          if opts.msd:
            m = MSD(natoms, dt, f1)

        f = Frame(natoms, n)
        f.parse_frame(buf)

        time.append(n*dt)
        if opts.stress:
          stress.append(f.stress)
          lat.append(f.lat)
        if opts.rdf:
          pairdist.update_rdf(f)
        if opts.vacf:
          c.update_vacf(n, f)
        if opts.msd:
          m.update_msd(n, f)
        continue
      if newframe:
        buf = ""
        newframe = False
      else:
        buf += line

  time = data['time']
  time = time - time[0]
  print()
  print("Analysing {} frames...".format(nframes))

# Plot the stress
  if opts.stress:
    stress = sp.array(stress)
    lat = sp.array(lat)
    mean_stress = sp.zeros((3,3))
    mean_lat = sp.zeros((3,3))
    for i in range(3):
      for j in range(3):
        mean_stress[i,j] = sp.mean(stress[:,i,j])
        mean_lat[i,j] = sp.mean(lat[:,i,j])
    plt.figure("Stress")

    if input_params['ensemble'][1] == "p":
      fig2, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
    else:
      fig2, (ax1,) = plt.subplots(nrows=1, ncols=1)

    plt.xlabel("t (fs)")
    ax1.set_ylabel("Stress (GPa)")
    ax2.set_ylabel("Cell dimension ($a_0$)")
    plt.xlim((time[opts.nskip], time[-1]))
    ax1.plot(time[opts.nskip:], stress[:,0,0], 'r-', label='xx', linewidth=1.0)
    ax1.plot(time[opts.nskip:], stress[:,1,1], 'g-', label='yy', linewidth=1.0)
    ax1.plot(time[opts.nskip:], stress[:,2,2], 'b-', label='zz', linewidth=1.0)
    ax1.plot((time[0],time[-1]), (mean_stress[0,0], mean_stress[0,0]), 'r-',
            label=r'$\langle S_{{xx}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,0]))
    ax1.plot((time[0],time[-1]), (mean_stress[1,1], mean_stress[1,1]), 'g-',
            label=r'$\langle S_{{yy}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,1]))
    ax1.plot((time[0],time[-1]), (mean_stress[2,2], mean_stress[2,2]), 'b-',
            label=r'$\langle S_{{zz}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,2]))

    if input_params['ensemble'][1] == "p":
      ax2.plot(time[opts.nskip:], lat[:,0,0], 'r-', label='a', linewidth=1.0)
      ax2.plot(time[opts.nskip:], lat[:,1,1], 'g-', label='b', linewidth=1.0)
      ax2.plot(time[opts.nskip:], lat[:,2,2], 'b-', label='c', linewidth=1.0)
      ax2.plot((time[0],time[-1]), (mean_lat[0,0], mean_lat[0,0]), 'r-',
              label=r'$\langle a \rangle$ = {0:<10.4f}'.format(mean_lat[0,0]))
      ax2.plot((time[0],time[-1]), (mean_lat[1,1], mean_lat[1,1]), 'g-',
              label=r'$\langle b \rangle$ = {0:<10.4f}'.format(mean_lat[1,1]))
      ax2.plot((time[0],time[-1]), (mean_lat[2,2], mean_lat[2,2]), 'b-',
              label=r'$\langle c \rangle$ = {0:<10.4f}'.format(mean_lat[2,2]))
      ax1.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
      ax2.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
      fig2.subplots_adjust(hspace=0)
      plt.setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)
      fig2.savefig("stress.pdf", bbox_inches='tight')

  # Plot the rdf
  if opts.rdf:
    pairdist.norm_rdf()
    pairdist.get_coordination()
    pairdist.plot_gr()
    if opts.dump:
      pairdist.dump_gr()

  # Plot the VACF
  if opts.vacf:
    c.norm_vacf()
    c.plot_vacf()
    if opts.dump:
      c.dump_vacf()

  # Plot the MSD
  if opts.msd:
    m.norm_msd()
    m.plot_msd()
    if opts.dump:
      m.dump_msd()
