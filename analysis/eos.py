#!/usr/local/bin/python3

import argparse
import scipy, re, sys, matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, leastsq
from pdb import set_trace
import warnings
warnings.filterwarnings("ignore")

class eos:

  def __init__(self, label, volume, energy):
    self.label = label
    self.volume = volume
    self.energy = energy

  def add_zpe(self, zpe):
    self.zpe = zpe
    self.helmholtz = self.energy + self.zpe

  def add_fine_grid(self, fine_grid):
    self.fine_grid = fine_grid

  def add_fit(self, e_fit, p_fit):
    self.e_fit = e_fit
    self.p_fit = p_fit

  def add_fit_zpe(self, f_fit, p_fit_zpe):
    self.f_fit = f_fit
    self.p_fit_zpe = p_fit_zpe

  def add_target(self, p, i):
    self.target_p = p
    self.target_i = i

  def add_target_zpe(self, p, i):
    self.target_p_zpe = p
    self.target_i_zpe = i

parser = argparse.ArgumentParser(description='Plot equation of state with and without zero-point energy correction', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-p', '--pressure', action='store', dest='press', default=None)
parser.add_argument('-i', '--infile', action='store', nargs='+', dest='infile', default='eos.dat', help='File containing volume/energy/ZPE data')
parser.add_argument('-V0', action='store', dest='V0', default=0.0)
parser.add_argument('-E0', action='store', dest='E0', default=0.0)
parser.add_argument('-B0', action='store', dest='B0', default=1.0)
parser.add_argument('-B0p', action='store', dest='B0p', default=6.0)
parser.add_argument('--birchm', action='store_true', dest='birchmfit', help='fit to Birch-Murnaghan equation of state')
parser.add_argument('--vinet', action='store_true', dest='vinetfit', help='fit to Vinet equation of state')
parser.add_argument('--quartic', action='store_true', dest='quartic', help='fit to quartic equation of state')
parser.add_argument('--zpe', action='store_true', dest='zpe', help='zero point energy column included in data file')
parser.add_argument('--au', action='store_true', dest='au', help='data in atomic units')
parser.add_argument('-v', '--volumes', action='store', nargs='+', dest='target_v', default=None, help='Compute E and P at these volumes')
cliopts = parser.parse_args()

if len(sys.argv) == 1:
  parser.print_help()
  sys.exit(0)

unit_vol = 1.0E-30
unit_ev = 1.60217653E-19
unit_giga = 1.0E9
eva32gpa = 160.21766208
habo32gpa = 29421.02648438959
ha2ev = 27.2113860217
bohr2ang = 0.529177249
const = unit_ev/unit_vol/unit_giga

tableheader = '{0: <16s}{1: <16s}{2: <16s}'
tablefmt = '{0: <16.6f}{1: <16.6f}{2:< 16.6f}'

npoints_fine = 100000
target_volume = None
if not cliopts.press == None:
  target_pressure = float(cliopts.press)
if not cliopts.target_v == None:
  target_volume = [float(v) for v in cliopts.target_v]

def eta(V, V0):
  return scipy.special.cbrt(V/V0)

def vinet(V, V0, E0, B0, B0p):
  e = eta(V, V0)
  f1 = -4*B0*V0/(B0p-1.)**2
  f2 = (3./2.)*(B0p-1.)*(1.-e)
  return E0 + f1*(1.-f2)*scipy.exp(f2)

def p_vinet(V, V0, E0, B0, B0p):
  e = eta(V, V0)
  f1 = (3./2.)*(B0p-1.)*(1.-e)
  return const*3.*B0*((1.-e)/e**2)*scipy.exp(f1)

def birchm(V, V0, E0, B0, B0p):
  e = eta(V, V0)
  f1 = (e**-2)-1
  f2 = 6-4*(e**-2)
  return E0 + (9*V0*B0/16)*((f1**3)*B0p + (f1**2)*f2)

def p_birchm(V, V0, E0, B0, B0p):
  e = eta(V, V0)
  return (3*B0p/2)*((e**-7)-(e**-5))*(1.0+0.75*(B0p-4.0)*((e**-2)-1))

def quartic(V, a0, b0, c0, d0, e0):
  return a0*V**4 + b0*V**3 + c0*V**2 + d0*V + e0

def p_quartic(V, a0, b0, c0, d0, e0):
  return -const*(4*a0*V**3 + 3*b0*V**2 + 2*c0*V + d0)

def quadraticfit(volumes, energies):
  a, b, c = scipy.polyfit(volumes, energies, 2)
  V0 = -b/(2*a)
  E0 = a*V0**2 + b*V0 + c
  B0 = 2*a*V0
  return V0, E0, B0

maxv = 0.0
minv = 1000000.0
eos_list = []
label = []
for inf in cliopts.infile:
  data = scipy.loadtxt(inf, dtype='float', skiprows=1)
  with open(inf, 'r') as infile:
    label = infile.readline().strip()
  volumes = data[:,0]
  energies = data[:,1]
  if cliopts.au:
    volumes = volumes*bohr2ang**3
    energies = energies*ha2ev
  vdiff = volumes[-1]-volumes[0]
  eos_i = eos(label, volumes, energies)
  if cliopts.zpe:
    zpes = data[:,2]
    if cliopts.au:
      zpes = zpes*ha2ev
    eos_i.add_zpe(zpes)
  if volumes[0] < minv:
    minv = volumes[0]
  if volumes[-1] > maxv:
    maxv = volumes[-1]
  v_min = volumes[0] - vdiff*0.1
  v_max = volumes[-1] + vdiff*0.1
  v_fine = scipy.linspace(volumes[0],volumes[-1], npoints_fine)
  eos_i.add_fine_grid(v_fine)
  eos_list.append(eos_i)

for eos_i in eos_list:

  if not cliopts.B0p:
    B0p = 4.
  else:
    B0p_init = float(cliopts.B0p)
  V0_init, E0_init, B0_init = quadraticfit(eos_i.volume, eos_i.energy)
  init_guess = (V0_init, E0_init, B0_init, B0p_init)
  if cliopts.vinetfit:

# Do a Vinet fit for the equation of state
    try:
      popt, pcov = curve_fit(vinet, eos_i.volume, eos_i.energy, p0=init_guess, maxfev=10000)
      if cliopts.zpe:
        popt2, pcov2 = curve_fit(vinet, eos_i.volume, eos_i.helmholtz, p0=popt, maxfev=10000)
    except RuntimeError:
      B0prange = scipy.linspace(2,8,num=25)
      for i in B0prange:
        popt2 = None
        init_guess = scipy.array([V0_init, E0_init, B0_init, i])
        try:
          popt, pcov = curve_fit(vinet, eos_i.volume, eos_i.energy, p0=init_guess, maxfev=10000)
          if cliopts.zpe:
            popt2, pcov2 = curve_fit(vinet, eos_i.volume, eos_i.helmholtz, p0=init_guess, maxfev=10000)
        except RuntimeError:
          continue
        if popt2.any():
          print("B0p = {}".format(i))
          print(popt2)
          break
      if not popt2.any():
        print("curve_fit failed")
        sys.exit(0)

    efit = vinet(v_fine, *popt)
    pfit_e = p_vinet(v_fine, *popt)
    eos_i.add_fit(efit, pfit_e)
    if cliopts.zpe:
      ffit = vinet(v_fine, *popt2)
      pfit_f = p_vinet(v_fine, *popt2)
      eos_i.add_fit_zpe(ffit, pfit_f)

  elif cliopts.quartic:
    c0, d0, e0 = scipy.polyfit(volumes, energies, 2)
    init_guess = scipy.array([0., 0., c0, d0, e0])
    popt, pcov = curve_fit(quartic, eos_i.volume, eos_i.energy, p0=init_guess, maxfev=10000)
    efit = quartic(v_fine, *popt)
    pfit_e = p_quartic(v_fine, *popt)
    eos_i.add_fit(efit, pfit_e)
    if cliopts.zpe:
      popt2, pcov2 = curve_fit(quartic, volumes, helmholtz, p0=popt, maxfev=10000)
      ffit = quartic(v_fine, *popt2)
      pfit_f = p_quartic(v_fine, *popt2)
      eos_i.add_fit_zpe(ffit, pfit_f)
  elif cliopts.birchmfit:
    popt, pcov = curve_fit(birchm, eos_i.volume, eos_i.energy, p0=init_guess, maxfev=10000)
    efit = birchm(v_fine, *popt)
    pfit_e = p_birchm(v_fine, *popt)
    eos_i.add_fit(efit, pfit_e)
    if cliopts.zpe:
      popt2, pcov2 = curve_fit(birchm, volumes, helmholtz, p0=popt, maxfev=10000)
      ffit = birchm(v_fine, *popt2)
      pfit_f = p_birchm(v_fine, *popt2)
      eos_i.add_fit_zpe(ffit, pfit_f)

  print(eos_i.label)
  print('Fitting parameters:')
  print('V0  = {0:<16.8f} +/- {1:<8.2g}'.format(popt[0], pcov[0,0]))
  print('E0  = {0:<16.8f} +/- {1:<8.2g}'.format(popt[1], pcov[1,1]))
  print('B0  = {0:<16.8f} +/- {1:<8.2g}'.format(popt[2], pcov[2,2]))
  print('B0p = {0:<16.8f} +/- {1:<8.2g}'.format(popt[3], pcov[3,3]))

  if target_volume:
    print('Energy/Pressure at specified volumes:')
    print(tableheader.format('V (A^3)', 'E (eV)', 'P (GPa)'))
    for v in target_volume:
      e = birchm(v, *popt)
      p = p_birchm(v, *popt)*eva32gpa
      print(tablefmt.format(v, e, p))
    print()

# Calculate position of line at target pressure
  if not cliopts.press == None:
    p_plot = 0.0
    p_target = target_pressure
    for i,v in enumerate(v_fine):
      diff = pfit_e[i]-p_plot
      if diff < p_target:
        p_target = pfit_e[i]
        i_target = i
        break
    try:
      pve = pfit_e[i_target] * v_fine[i_target] / eva32gpa
    except:
      set_trace()
    eos_i.add_target(p_target, i_target)

    if cliopts.zpe:
      p_target_f = target_pressure
      for i,v in enumerate(v_fine):
        diff = pfit_f[i]-p_plot
        if diff < p_target_f:
          p_target_f = pfit_f[i]
          i_target_f = i
          break
      pvf = pfit_f[i_target_f] * v_fine[i_target_f] / eva32gpa
      eos_i.add_target_zpe(p_target, i_target)


    print('Without vibrations: P = {0:>6.2f}, V = {1:>12.8f}, PV = {2:>10.6f}'.format(pfit_e[i_target], v_fine[i_target], pve))
    print('                    E = {0:>10.6f}, H = {1:>10.6f}'.format(efit[i_target], efit[i_target]+pve))
    if cliopts.zpe:
      print('With vibrations:    P = {0:>6.2f}, V = {1:>12.8f}, PV = {2:>10.6f}'.format(pfit_f[i_target_f], v_fine[i_target_f], pvf))
      print('                    E = {0:>10.6f}, H = {1:>10.6f}'.format(ffit[i_target_f], ffit[i_target_f]+pvf))

  e_min = energies[-1] - scipy.absolute((energies[0]-energies[-1])*1.0)
  e_max = energies[0] + scipy.absolute((energies[0]-energies[-1])*1.0)

# Plot graphs
plt.ioff()
ax1 = plt.subplot(211)
#  l = 'E = {0:>6.2f}, V = {1:>8.4f}'.format(efit[i_target],
#                                            v_fine[i_target])
plt.ylabel('E (eV)')
plt.setp(ax1.get_xticklabels(), visible=False)
#plt.ylim(e_min,e_max)
for e in eos_list:
  # Plot the data points
  plt.plot(e.volume, e.energy, 'ko', markerfacecolor='None')
  # Plot the fit
  plt.plot(e.fine_grid, e.e_fit, linestyle='-', label=e.label)
  if cliopts.zpe:
    plt.plot(e.volume, e.helmholtz, 'ko', markerfacecolor='None', markeredgewidth=1.5)
    plt.plot(e.fine_grid, e.f_fit, linestyle='-', label=e.label)

ax1.autoscale(False)
plt.gca().set_color_cycle(None)
for e in eos_list:
  if not cliopts.press == None:
    # plot a line for the target pressure
    plt.plot((e.fine_grid[e.target_i],e.fine_grid[e.target_i]), (-100000,1000000), linestyle='--')
plt.legend(loc='lower left')

ax2 = plt.subplot(212, sharex=ax1)
#l1 = 'Static:     P = {0:>6.2f}, V = {1:>8.4f}'.format(pfit_e[i_target], v_fine[i_target])
plt.xlabel(r'V ($\AA ^3$)')
plt.ylabel('P (GPa)')
for e in eos_list:
  plt.plot(e.fine_grid, e.p_fit*eva32gpa, linestyle='-')
  if cliopts.zpe:
#    l2 = 'Vibrations: P = {0:>6.2f}, V = {1:>8.4f}'.format(pfit_f[i_target_f], v_fine[i_target_f])
    plt.plot(e.fine_grid, e.p_fit_zpe, linestyle='-')

ax2.autoscale(False)
plt.gca().set_color_cycle(None)
for e in eos_list:
  plt.plot((e.fine_grid[e.target_i],e.fine_grid[e.target_i]), (-100000,1000000), linestyle='--')

if not cliopts.press == None:
  ax2.autoscale(False)
  plt.plot((minv, maxv), (target_pressure, target_pressure), 'k', linestyle = '--')

plt.legend(loc='lower left')
plt.xlim(minv,maxv)
#plt.ylim(ymin=0)
plt.tight_layout()
plt.savefig('eos.pdf', bbox_inches='tight')
