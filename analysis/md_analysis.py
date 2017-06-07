#!/usr/local/bin/python3

import argparse
import re
import scipy as sp
import matplotlib.pyplot as plt
from frame import Frame


frame_re = re.compile('frame')
endframe_re = re.compile('end frame')
cell_re = re.compile('cell_vectors(.*?)end cell_vectors', re.M | re.S)
stress_re = re.compile('stress_tensor(.*?)end stress_tensor', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
position_re = re.compile('positions(.*?)end positions', re.M | re.S)
velocity_re = re.compile('velocities(.*?)end velocities', re.M | re.S)
force_re = re.compile('forces(.*?)end forces', re.M | re.S)


def parse_frame(buf, f):
  m = re.search(cell_re, buf)
  lines = m.group(1).strip().splitlines()
  for i in range(3):
    bits = lines[i].strip().split()
    for j in range(3):
      f.lat[i,j] = float(bits[j])

  m = re.search(stress_re, buf)
  lines = m.group(1).strip().splitlines()
  for i in range(3):
    bits = lines[i].strip().split()
    for j in range(3):
      f.stress[i,j] = float(bits[j])

  m = re.search(position_re, buf)
  lines = m.group(1).strip().splitlines()
  nat = len(lines)
  for i in range(nat):
    bits = lines[i].strip().split()
    bits.pop(0)
    f.species[i] = int(bits.pop(0))
    for j in range(3):
      f.r[i,j] = float(bits[j])

  m = re.search(velocity_re, buf)
  lines = m.group(1).strip().splitlines()
  nat = len(lines)
  for i in range(nat):
    bits = lines[i].strip().split()
    bits.pop(0)
    bits.pop(0)
    for j in range(3):
      f.v[i,j] = float(bits[j])

  m = re.search(force_re, buf)
  lines = m.group(1).strip().splitlines()
  nat = len(lines)
  for i in range(nat):
    bits = lines[i].strip().split()
    bits.pop(0)
    bits.pop(0)
    for j in range(3):
      f.f[i,j] = float(bits[j])


parser = argparse.ArgumentParser(description='Analyse a MD trajectory',
                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-d', '--dump', action='store', dest='dumpfile',
                    default='dump.out', help='MD dump file')
parser.add_argument('-s', '--stat', action='store', dest='statfile',
                    default='stat.out', help='MD statistics file')
parser.add_argument('-m', '--mdin', action='store', dest='mdinfile',
                    default='md.in', help='MD parameters file')
parser.add_argument('--skip', action='store', dest='nskip', default=0,
                    type=int, help='Number of equilibration steps to skip')

opts = parser.parse_args()

# Parse the md.in parameters file
mdin_params = {}
with open(opts.mdinfile, 'r') as infile:
  for line in infile:
    param, value = line.split()
    mdin_params[param] = value

natoms = int(mdin_params['natoms'])
dt = float(mdin_params['dt'])

# Parse the stat.out statistics file
step = []
pe = []
ke = []
nhc = []
E = []
T = []
P = []
with open(opts.statfile, 'r') as statfile:
  statfile.readline()
  for line in statfile:
    if mdin_params['ensemble'] == 'nve':
      a, b, c, d, e, f = line.strip().split()
      step.append(int(a))
      pe.append(float(b))
      ke.append(float(c))
      E.append(float(d))
      T.append(float(e))
      P.append(float(f))
    elif mdin_params['thermo_type'] == 'nhc':
      a, b, c, d, e, f, g = line.strip().split()
      step.append(int(a))
      pe.append(float(b))
      ke.append(float(c))
      nhc.append(float(d))
      E.append(float(e))
      T.append(float(f))
      P.append(float(g))

step = sp.array(step)
pe = sp.array(pe)
ke = sp.array(ke)
nhc = sp.array(nhc)
E = sp.array(E)
T = sp.array(T)
P = sp.array(P)
E_avg = sp.mean(E[opts.nskip:-1])
T_avg = sp.mean(T[opts.nskip:-1])
P_avg = sp.mean(P[opts.nskip:-1])

# Plot the statistics
fig1, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=True)

plt.xlim((opts.nskip,step[-1]))
ax1.plot(step, pe, 'r-', label='Potential energy')
ax1.plot(step, ke, 'b-', label='Kinetic energy')
if mdin_params['thermo_type'] == 'nhc':
  ax1.plot(step, nhc, 'g-', label='NHC energy')
ax2.plot(step, E)
ax2.plot((opts.nskip,step[-1]), (E_avg,E_avg), '-',
      label=r'$\langle E \rangle$ = {0:<12.4f}'.format(E_avg))
ax3.plot(step, T)
ax3.plot((opts.nskip,step[-1]), (T_avg,T_avg), '-',
      label=r'$\langle T \rangle$ = {0:<12.4f}'.format(T_avg))
ax4.plot(step, P)
ax4.plot((opts.nskip,step[-1]), (P_avg,P_avg), '-',
      label=r'$\langle P \rangle$ = {0:<12.4f}'.format(P_avg))
ax1.set_ylabel("E")
ax2.set_ylabel("E")
ax3.set_ylabel("T")
ax4.set_ylabel("P")
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()

plt.xlabel("MD step")
fig1.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)
fig1.savefig("stats.pdf", bbox_inches='tight')

# Parse the dump.out dump file
nframes = 0
newframe = True
nbins = 100
buf = ""
time = []
vacf = []
msd = []
stress = []
sdistr = sp.zeros(nbins, dtype='float')
vdistr = sp.zeros((nbins,3), dtype='float')
with open(opts.dumpfile, 'r') as dumpfile:
  while True:
    line = dumpfile.readline()
    if not line:
      break
    if re.match(frame_re, line):
      n = int(line.split()[1])

    if re.match(endframe_re, line):
      newframe = True
      if n < opts.nskip:
        continue
      else:
        nframes += 1
      if nframes == 1:
        f1 = Frame(natoms,n)
        parse_frame(buf,f1)
      f = Frame(natoms, n)
      parse_frame(buf, f)
      stress.append(f.stress)
      time.append(n*dt)
      vacf.append(f.update_vacf(f1))
      msd.append(f.update_msd(f1))
      sdistr_tmp, bin_edges = f.update_sdistr(nbins)
      vdistr_tmp, bin_edges = f.update_vdistr(nbins)
      sdistr += sdistr_tmp
      vdistr += vdistr_tmp
      continue
    if newframe:
      buf = ""
      newframe = False
    else:
      buf += line

vacf = sp.array(vacf)
msd = sp.array(msd)
stress = sp.array(stress)
time = sp.array(time)
time = time - time[0]
vdistr = vdistr/float(nframes)
sdistr = sdistr/float(nframes)

# Plot the stress
mean_stress = sp.zeros((3,3))
for i in range(3):
  for j in range(3):
    mean_stress[i,j] = sp.mean(stress[:,i,j])
plt.figure("Stress")
fig2, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, sharex=True)
plt.xlabel("t")
ax1.set_ylabel("Stress")
ax2.set_ylabel("Stress")
plt.xlim((time[0], time[-1]))
ax1.plot(time, stress[:,0,0], 'r-', label='xx')
ax1.plot(time, stress[:,1,1], 'g-', label='yy')
ax1.plot(time, stress[:,2,2], 'b-', label='zz')
ax1.plot((time[0],time[-1]), (mean_stress[0,0], mean_stress[0,0]), 'r-',
        label=r'$\langle S_{{xx}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,0]))
ax1.plot((time[0],time[-1]), (mean_stress[1,1], mean_stress[1,1]), 'g-',
        label=r'$\langle S_{{yy}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,1]))
ax1.plot((time[0],time[-1]), (mean_stress[2,2], mean_stress[2,2]), 'b-',
        label=r'$\langle S_{{zz}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,2]))

ax2.plot(time, stress[:,0,1], 'r-', label='xy')
ax2.plot(time, stress[:,1,0], 'r--', label='yx')
ax2.plot(time, stress[:,1,2], 'g-', label='yz')
ax2.plot(time, stress[:,2,1], 'g--', label='zy')
ax2.plot(time, stress[:,0,2], 'b-', label='xz')
ax2.plot(time, stress[:,2,0], 'b--', label='zx')
ax2.plot((time[0],time[-1]), (mean_stress[0,1], mean_stress[0,1]), 'r-',
        label=r'$\langle S_{{xy}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,1]))
ax2.plot((time[0],time[-1]), (mean_stress[1,0], mean_stress[1,0]), 'r--',
        label=r'$\langle S_{{yx}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,0]))
ax2.plot((time[0],time[-1]), (mean_stress[0,2], mean_stress[0,2]), 'g-',
        label=r'$\langle S_{{xz}} \rangle$ = {0:<10.4f}'.format(mean_stress[0,2]))
ax2.plot((time[0],time[-1]), (mean_stress[2,0], mean_stress[2,0]), 'g--',
        label=r'$\langle S_{{zx}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,0]))
ax2.plot((time[0],time[-1]), (mean_stress[1,2], mean_stress[1,2]), 'b-',
        label=r'$\langle S_{{yz}} \rangle$ = {0:<10.4f}'.format(mean_stress[1,2]))
ax2.plot((time[0],time[-1]), (mean_stress[2,1], mean_stress[2,1]), 'b--',
        label=r'$\langle S_{{zy}} \rangle$ = {0:<10.4f}'.format(mean_stress[2,1]))
ax1.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
ax2.legend(bbox_to_anchor=(1.05,1), loc=2, borderaxespad=0.)
fig2.subplots_adjust(hspace=0)
plt.setp([a.get_xticklabels() for a in fig1.axes[:-1]], visible=False)
fig2.savefig("stress.pdf", bbox_inches='tight')


# Plot the VACF
plt.figure("VACF")
plt.xlabel("t")
plt.ylabel("C(t)")
plt.xlim((time[0],time[-1]))
plt.plot(time, vacf)
plt.plot((0,time[-1]), (0, 0), 'k-')
plt.savefig("vacf.pdf", bbox_inches='tight')

# Plot the MSD
plt.figure("MSD")
plt.xlabel("t")
plt.ylabel("MSD(t)")
plt.xlim((time[0],time[-1]))
plt.plot(time, msd)
plt.plot((0,time[-1]), (0, 0), 'k-')
plt.ylim(ymin=0)
plt.savefig("msd.pdf", bbox_inches='tight')

# Plot the velocity distribution
plt.figure("v_distribution")
fig3, ((ax1, ax2), (ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
ax1.hist(sdistr, bins=bin_edges, label='speed')
ax2.hist(vdistr[:,0], bins=bin_edges, label=r'$v_x$')
ax3.hist(vdistr[:,1], bins=bin_edges, label=r'$v_y$')
ax4.hist(vdistr[:,2], bins=bin_edges, label=r'$v_z$')
plt.savefig("vdistr.pdf", bbox_inches='tight')
