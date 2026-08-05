#!/usr/bin/env python3
"""
plot_fast_raw.py — Full visualization suite for FAST raw FRO HDF5 files.

Produces all scientifically meaningful plots from a FAST raw .fro file:
  1. Mean spectrum full band (log scale, XX+YY overlaid)
  2. Mean spectrum zoomed on HI line (XX+YY overlaid, with scientific notes)
  3. Waterfall (time vs velocity, HI band, mean of XX+YY)
  4. Coverage map (RA/Dec of all pointings, coloured by elevation)
  5. Moment-0 map (interpolated on spatial grid, mean XX+YY)
  6. Position-velocity diagram along RA axis (mean XX+YY)
  7. Elevation vs time
  8. XX vs YY polarization comparison (HI band)

Note: data are raw uncalibrated power (counts). No baseline subtraction
or flux calibration applied. Spectral features reflect receiver bandpass
and RFI, not calibrated brightness temperature.

Author: Francesco Di Giovanni (Claude AI assisted) - FRO/HALO Project
"""

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import h5py

C_MS = 299792458.0
HI_REST_HZ = 1420405752.0


def load_fro(fro_path):
    with h5py.File(fro_path, 'r') as f:
        freq_hz  = f['Spectra/frequency_axis_hz'][:]
        spectra  = f['Spectra/raw_spectrum'][:]
        pol      = f['Spectra/polarization'][:].astype(str)
        ra       = f['Pointing/ra_deg'][:]
        dec      = f['Pointing/dec_deg'][:]
        el       = f['Pointing/elevation_deg'][:]
        unix     = f['Time/unix_time'][:]
    return dict(freq_hz=freq_hz, spectra=spectra, pol=pol,
                ra=ra, dec=dec, el=el, unix=unix)


def freq_to_vel(freq_hz):
    return (HI_REST_HZ - freq_hz) / HI_REST_HZ * C_MS / 1000.0


def hi_mask(freq_hz):
    return (freq_hz >= 1380e6) & (freq_hz <= 1430e6)


def mean_both_pol(d, mask_freq=None):
    """Return mean spectrum averaged over both polarizations."""
    xx = d['pol'] == 'XX'
    yy = d['pol'] == 'YY'
    if mask_freq is not None:
        sxx = np.nanmean(d['spectra'][xx][:, mask_freq], axis=0)
        syy = np.nanmean(d['spectra'][yy][:, mask_freq], axis=0)
    else:
        sxx = np.nanmean(d['spectra'][xx], axis=0)
        syy = np.nanmean(d['spectra'][yy], axis=0)
    return sxx, syy


def style_ax(ax):
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_color('white')


def plot_mean_spectrum_full(d, target, plots_dir):
    freq = d['freq_hz']
    sxx, syy = mean_both_pol(d)
    fig, ax = plt.subplots(figsize=(12, 5))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.semilogy(freq / 1e6, sxx,  color='cyan',   linewidth=0.5, label='XX', alpha=0.9)
    ax.semilogy(freq / 1e6, syy,  color='orange', linewidth=0.5, label='YY', alpha=0.9)
    ax.axvline(HI_REST_HZ / 1e6, color='yellow', linestyle='--',
               linewidth=0.8, label='HI rest 1420.4 MHz')
    ax.set_xlabel('Frequency (MHz)', color='white')
    ax.set_ylabel('Power (counts, log scale)', color='white')
    ax.set_title(f'{target} — FAST Mean Spectrum full band (XX+YY, log scale)', color='white')
    ax.legend(facecolor='black', labelcolor='white')
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_spectrum_full.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def plot_hi_spectrum(d, target, plots_dir):
    freq = d['freq_hz']
    mask = hi_mask(freq)
    vel  = freq_to_vel(freq[mask])
    sxx, syy = mean_both_pol(d, mask_freq=mask)
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.plot(vel, sxx, color='cyan',   linewidth=1.0, label='XX', alpha=0.9)
    ax.plot(vel, syy, color='orange', linewidth=1.0, label='YY', alpha=0.9)
    ax.axvline(-180, color='yellow', linestyle='--', linewidth=0.8,
               label='M33 v_sys ~-180 km/s')
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.8,
               label='Galactic HI ~0 km/s (local)')
    ax.set_xlabel('LSR velocity (km/s)', color='white')
    ax.set_ylabel('Power (counts, uncalibrated)', color='white')
    ax.set_title(f'{target} — FAST HI band spectrum (XX+YY)\n'
                 f'Note: raw data — peak near 0 km/s is local Galactic HI; '
                 f'M33 at -180 km/s requires calibration+baseline subtraction',
                 color='white', fontsize=9)
    ax.legend(facecolor='black', labelcolor='white')
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_spectrum_hi.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def plot_waterfall(d, target, plots_dir):
    freq = d['freq_hz']
    mask = hi_mask(freq)
    vel  = freq_to_vel(freq[mask])
    xx = d['pol'] == 'XX'
    yy = d['pol'] == 'YY'
    # Mean of XX and YY per integration (zip by index)
    sxx = d['spectra'][xx][:, mask]
    syy = d['spectra'][yy][:, mask]
    n = min(sxx.shape[0], syy.shape[0])
    wf = (sxx[:n] + syy[:n]) / 2.0
    fig, ax = plt.subplots(figsize=(10, 7))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    vmin = np.nanpercentile(wf, 2)
    vmax = np.nanpercentile(wf, 98)
    im = ax.imshow(wf, origin='upper', aspect='auto', cmap='inferno',
                   vmin=vmin, vmax=vmax,
                   extent=[vel[0], vel[-1], n, 0])
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('Power (counts)', color='white')
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax.axvline(-180, color='yellow', linestyle='--', linewidth=0.8,
               label='M33 v_sys')
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.8,
               label='Galactic HI')
    ax.set_xlabel('LSR velocity (km/s)', color='white')
    ax.set_ylabel('Integration index', color='white')
    ax.set_title(f'{target} — FAST HI Waterfall (mean XX+YY)', color='white')
    ax.legend(facecolor='black', labelcolor='white')
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_waterfall.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def plot_coverage(d, target, plots_dir):
    xx = d['pol'] == 'XX'
    fig, ax = plt.subplots(figsize=(7, 6))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    sc = ax.scatter(d['ra'][xx], d['dec'][xx],
                    c=d['el'][xx], cmap='plasma', s=4, alpha=0.7)
    cb = plt.colorbar(sc, ax=ax)
    cb.set_label('Elevation (deg)', color='white')
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax.set_xlabel('RA (deg)', color='white')
    ax.set_ylabel('Dec (deg)', color='white')
    ax.set_title(f'{target} — FAST OTF Coverage Map', color='white')
    ax.invert_xaxis()
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_map_coverage.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def plot_pv_diagram(d, target, plots_dir):
    freq = d['freq_hz']
    mask = hi_mask(freq)
    vel  = freq_to_vel(freq[mask])
    xx = d['pol'] == 'XX'
    yy = d['pol'] == 'YY'
    sxx = d['spectra'][xx][:, mask]
    syy = d['spectra'][yy][:, mask]
    n   = min(sxx.shape[0], syy.shape[0])
    mean_spec = (sxx[:n] + syy[:n]) / 2.0
    ra  = d['ra'][xx][:n]
    sort_idx     = np.argsort(ra)
    specs_sorted = mean_spec[sort_idx]
    ra_offset    = (ra[sort_idx] - ra[sort_idx].mean()) * 60.0
    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    vmin = np.nanpercentile(specs_sorted, 2)
    vmax = np.nanpercentile(specs_sorted, 98)
    im = ax.imshow(specs_sorted.T, origin='lower', aspect='auto',
                   cmap='inferno', vmin=vmin, vmax=vmax,
                   extent=[ra_offset[0], ra_offset[-1], vel[0], vel[-1]])
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('Power (counts)', color='white')
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax.axhline(-180, color='yellow', linestyle='--', linewidth=0.8,
               label='M33 v_sys ~-180 km/s')
    ax.axhline(0, color='gray', linestyle=':', linewidth=0.8,
               label='Galactic HI ~0 km/s')
    ax.set_xlabel('RA offset (arcmin)', color='white')
    ax.set_ylabel('LSR velocity (km/s)', color='white')
    ax.set_title(f'{target} — FAST HI Position-Velocity (mean XX+YY)', color='white')
    ax.legend(facecolor='black', labelcolor='white')
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_pv_diagram.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def plot_elevation_time(d, target, plots_dir):
    xx   = d['pol'] == 'XX'
    unix = d['unix'][xx]
    el   = d['el'][xx]
    t_min = (unix - unix[0]) / 60.0
    fig, ax = plt.subplots(figsize=(10, 4))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.plot(t_min, el, color='lime', linewidth=1.0)
    ax.set_xlabel('Time from start (min)', color='white')
    ax.set_ylabel('Elevation (deg)', color='white')
    ax.set_title(f'{target} — FAST Elevation vs Time', color='white')
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_elevation_time.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def plot_xx_yy_comparison(d, target, plots_dir):
    freq = d['freq_hz']
    mask = hi_mask(freq)
    vel  = freq_to_vel(freq[mask])
    sxx, syy = mean_both_pol(d, mask_freq=mask)
    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.plot(vel, sxx, color='cyan',   linewidth=1.0, label='XX', alpha=0.9)
    ax.plot(vel, syy, color='orange', linewidth=1.0, label='YY', alpha=0.9)
    ax.axvline(-180, color='yellow', linestyle='--', linewidth=0.8,
               label='M33 v_sys ~-180 km/s')
    ax.axvline(0, color='gray', linestyle=':', linewidth=0.8,
               label='Galactic HI ~0 km/s')
    ax.set_xlabel('LSR velocity (km/s)', color='white')
    ax.set_ylabel('Power (counts)', color='white')
    ax.set_title(f'{target} — FAST HI XX vs YY Polarization Comparison', color='white')
    ax.legend(facecolor='black', labelcolor='white')
    style_ax(ax)
    plt.tight_layout()
    out = os.path.join(plots_dir, f'{target}_spectrum_xx_yy.png')
    plt.savefig(out, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {out}")


def main():
    if len(sys.argv) > 1:
        fro_path = sys.argv[1]
    else:
        fro_path = input("Path to FAST .fro file: ").strip()
    fro_path = os.path.expanduser(fro_path)
    if not os.path.isfile(fro_path):
        print(f"ERROR: file not found: {fro_path}")
        sys.exit(1)

    print(f"Loading {fro_path} ...")
    d = load_fro(fro_path)
    surveys_root = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Surveys')
    target    = 'M33'
    plots_dir = os.path.join(surveys_root, target, 'Plots')
    os.makedirs(plots_dir, exist_ok=True)

    print(f"Target   : {target}")
    print(f"Spectra  : {d['spectra'].shape[0]}")
    print(f"Channels : {d['spectra'].shape[1]}")
    print(f"Plots dir: {plots_dir}")

    print("\n1. Mean spectrum full band (log scale) ...")
    plot_mean_spectrum_full(d, target, plots_dir)

    print("2. HI spectrum zoomed ...")
    plot_hi_spectrum(d, target, plots_dir)

    print("3. Waterfall (HI band) ...")
    plot_waterfall(d, target, plots_dir)

    print("4. Coverage map ...")
    plot_coverage(d, target, plots_dir)


    print("6. PV diagram ...")
    plot_pv_diagram(d, target, plots_dir)

    print("7. Elevation vs time ...")
    plot_elevation_time(d, target, plots_dir)

    print("8. XX vs YY comparison ...")
    plot_xx_yy_comparison(d, target, plots_dir)

    print(f"\nAll plots saved in: {plots_dir}")
    for f in sorted(os.listdir(plots_dir)):
        if f.endswith('.png'):
            subprocess.Popen(['xdg-open', os.path.join(plots_dir, f)])


if __name__ == '__main__':
    main()