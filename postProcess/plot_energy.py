#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import argparse
import os

plt.rcParams['text.usetex'] = True
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Computer Modern']
plt.rcParams['font.size'] = 12
plt.rcParams['axes.labelsize'] = 14
plt.rcParams['axes.titlesize'] = 16
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 12
plt.rcParams['ytick.labelsize'] = 12
plt.rcParams['figure.dpi'] = 150
plt.rcParams['savefig.dpi'] = 300
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['savefig.pad_inches'] = 0.1

def read_energy_data(filename):
  """Read energy data from the simulation output file."""
  data = np.loadtxt(filename, skiprows=1)
  
  t = data[:, 1]
  ke = data[:, 2]
  se = data[:, 3]
  _ = data[:, 4] # dissipation rate -- do not plot!
  Ediss = data[:, 5]
  Etot = data[:, 6]
  
  return t, ke, se, Ediss, Etot

def plot_energy(filename):
  """Create energy plot with nice aesthetics."""
  t, ke, se, Ediss, Etot = read_energy_data(filename)
  
  fig, ax = plt.subplots(figsize=(10, 6))
  
  ax.plot(t, ke, 'b-', linewidth=2, label=r'$E_{\mathrm{kinetic}}$')
  ax.plot(t, se, 'r-', linewidth=2, label=r'$E_{\mathrm{surface}}$')
  ax.plot(t, Ediss, 'g-', linewidth=2, label=r'$E_{\mathrm{dissipated}}$')
  ax.plot(t, Etot, 'k--', linewidth=2.5, label=r'$E_{\mathrm{total}}$')
  
  ax.set_xlabel(r'Time $t$', fontsize=14)
  ax.set_ylabel(r'Energy', fontsize=14)
  ax.set_title(r'Energy Evolution During Drop Impact', fontsize=16)
  
  ax.grid(True, alpha=0.3, linestyle='--')
  ax.legend(loc='best', frameon=True, fancybox=True, shadow=True)
  
  # Keep all spines visible for box frame
  for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(1.5)
  ax.tick_params(direction='out', length=6, width=1)

  ax.set_xlim(0, t[-1])
  ax.set_ylim(0, 2.2)
  
  plt.tight_layout()
  
  output_file = os.path.splitext(filename)[0] + '_plot.png'
  plt.savefig(output_file, bbox_inches='tight')
  print(f"Plot saved to: {output_file}")
  
  plt.show()

def main():
  parser = argparse.ArgumentParser(description='Plot energy data from Basilisk simulation')
  parser.add_argument('filename', nargs='?', 
                      default='../simulationCases/3-DropImpactOnSolids/energy.dat',
                      help='Input energy data file (default: ../simulationCases/3-DropImpactOnSolids/energy.dat)')
  
  args = parser.parse_args()
  
  if not os.path.exists(args.filename):
    print(f"Error: File '{args.filename}' not found.")
    return
  
  plot_energy(args.filename)

if __name__ == '__main__':
  main()