#!/usr/bin/env python3
# plot_gbt.py
# Plotting script for GBT SDFITS .fro files (FRO HDF5 v17)
# Produces: averaged spectrum per polarization + all spectra overlay
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy
# 2026-07-29

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py

RESTFREQ_HZ = 1420.405750e6  # HI rest frequency


def vel_axis(freq_hz):
    return (RESTFREQ_HZ - freq_hz) / RESTFREQ_HZ * 2.998e5


def load_fro(path):
    import re
    with h5py.File(path, 'r') as f:
        freq_hz = f['Spectra/frequency_axis_hz'][:]
        spectra = f['Spectra/raw_spectrum'][:]
        pols = f['Spectra/polarization'][:].astype(str)
        target = f['Observation'].attrs.get('target',
                 f['Source'].attrs.get('target', os.path.basename(path)))
        facility = f['Source'].attrs.get('facility', 'GBT')
        obs_id = f['Observation'].attrs.get('observation_id', '')
        # Read RESTFREQ from Header notes
        notes = f['Header'].attrs.get('notes', '')
        m = re.search(r'RESTFREQ\s+([\d.]+)\s+Hz', str(notes))
        restfreq_hz = float(m.group(1)) if m else 1420.405750e6
    return freq_hz, spectra, pols, str(target), str(facility), str(obs_id), restfreq_hz


def plot_fro(fro_path, out_dir):
    freq_hz, spectra, pols, target, facility, obs_id, restfreq_hz = load_fro(fro_path)
    vel = (restfreq_hz - freq_hz) / restfreq_hz * 2.998e5
    tag = os.path.splitext(os.path.basename(fro_path))[0]
    unique_pols = list(dict.fromkeys(pols))
    colors = {'XX': '#1f77b4', 'YY': '#d62728',
              'RR': '#1f77b4', 'LL': '#d62728', 'I': '#2ca02c'}
    title = f'GBT — {target}' if target and target != '?' else f'GBT — {tag}'

    # Figure 1: averaged spectrum per polarization
    fig, ax = plt.subplots(figsize=(11, 5), facecolor='black')
    ax.set_facecolor('black')
    for pol in unique_pols:
        mask = pols == pol
        avg = np.nanmean(spectra[mask], axis=0)
        ax.plot(vel, avg, color=colors.get(pol, '#7f7f7f'),
                linewidth=0.8, label=pol)
    ax.axvline(0, color='yellow', linewidth=0.7, linestyle='--',
               alpha=0.6, label='v=0 km/s')
    ax.set_xlabel('LSR Velocity (km/s)', color='white', fontsize=11)
    ax.set_ylabel('T$_A$ (K)', color='white', fontsize=11)
    ax.set_title(title, color='white', fontsize=13, pad=10)
    ax.tick_params(colors='white')
    ax.spines[:].set_edgecolor('white')
    ax.legend(fontsize=9, facecolor='#222', edgecolor='white', labelcolor='white')
    ax.grid(True, color='#333', linewidth=0.4)
    fig.tight_layout()
    out1 = os.path.join(out_dir, f'{tag}_avg_spectrum.png')
    fig.savefig(out1, dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    print(f'  Written: {out1}')

    # Figure 2: all spectra overlay
    fig, ax = plt.subplots(figsize=(11, 5), facecolor='black')
    ax.set_facecolor('black')
    for spec, pol in zip(spectra, pols):
        ax.plot(vel, spec, color=colors.get(pol, '#7f7f7f'),
                linewidth=0.4, alpha=0.5)
    ax.axvline(0, color='yellow', linewidth=0.7, linestyle='--',
               alpha=0.6, label='v=0 km/s')
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=colors.get(p, '#7f7f7f'),
                      linewidth=1.5, label=p) for p in unique_pols]
    handles.append(Line2D([0], [0], color='yellow', linewidth=1,
                          linestyle='--', label='v=0 km/s'))
    ax.legend(handles=handles, fontsize=9, facecolor='#222',
              edgecolor='white', labelcolor='white')
    ax.set_xlabel('LSR Velocity (km/s)', color='white', fontsize=11)
    ax.set_ylabel('T$_A$ (K)', color='white', fontsize=11)
    ax.set_title(f'{title} — all spectra', color='white', fontsize=13, pad=10)
    ax.tick_params(colors='white')
    ax.spines[:].set_edgecolor('white')
    ax.grid(True, color='#333', linewidth=0.4)
    fig.tight_layout()
    out2 = os.path.join(out_dir, f'{tag}_all_spectra.png')
    fig.savefig(out2, dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    print(f'  Written: {out2}')
    return [out1, out2]


def main():
    if len(sys.argv) < 2:
        base = os.path.expanduser('~/FRO/Parsers/GBT/Surveys/')
        fro_files = sorted(glob.glob(
            os.path.join(base, '*', 'GBT_FRO_*', '*.h5')))
        if not fro_files:
            print('No .fro files found. Usage: python3 plot_gbt.py <file.h5> [out_dir]')
            sys.exit(1)
        print(f'Found {len(fro_files)} .fro file(s):')
        for i, f in enumerate(fro_files):
            print(f'  [{i}] {os.path.basename(f)}')
        choice = input('Select file number (or Enter for all): ').strip()
        selected = fro_files if choice == '' else [fro_files[int(choice)]]
    else:
        selected = [os.path.expanduser(sys.argv[1])]

    import subprocess
    all_plots = []
    for fro_path in selected:
        fro_dir = os.path.dirname(fro_path)
        obs_dir = os.path.dirname(fro_dir)
        plots_dirname = os.path.basename(fro_dir).replace('GBT_FRO_', 'GBT_Plots_')
        out_dir = os.path.expanduser(sys.argv[2]) if len(sys.argv) >= 3 else \
                  os.path.join(obs_dir, plots_dirname)
        os.makedirs(out_dir, exist_ok=True)
        print(f'\nPlotting: {os.path.basename(fro_path)}')
        all_plots.extend(plot_fro(fro_path, out_dir))
    for p in all_plots:
        subprocess.Popen(['xdg-open', p], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print('\nDone.')


if __name__ == '__main__':
    main()
