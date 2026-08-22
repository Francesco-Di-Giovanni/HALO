#!/usr/bin/env python3
# seti_plot_salsa.py
# HALO/FRO Project - SETI narrowband analysis for SALSA Raw FITS
# Author: Francesco Di Giovanni (Claude AI assisted) - FRO/HALO Project
# Usage: python3 seti_plot_salsa.py <path_to_fits>
#
# Method:
#   1. Bandpass removal via running median (window=401 ch, ~49 kHz)
#   2. Local noise estimation via rolling MAD (window=501 ch, ~61 kHz)
#   3. Narrowband peak search: |SNR| > threshold AND width <= 5 channels
#   4. Three diagnostic plots saved to ../Plots/ relative to FITS location

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import median_filter, generic_filter

# ── input ──────────────────────────────────────────────────────────────────────
if len(sys.argv) < 2:
    print("Usage: python3 seti_plot_salsa.py <path_to_fits>")
    sys.exit(1)

fits_path = sys.argv[1]
if not os.path.isfile(fits_path):
    print(f"ERROR: file not found: {fits_path}")
    sys.exit(1)

# Output dir: ../Plots/ relative to FITS file
fits_dir  = os.path.dirname(os.path.abspath(fits_path))
epoch_dir = os.path.dirname(fits_dir)
plots_dir = os.path.join(epoch_dir, 'Plots')
os.makedirs(plots_dir, exist_ok=True)
print(f"FITS  : {fits_path}")
print(f"Plots : {plots_dir}")

# ── load FITS ──────────────────────────────────────────────────────────────────
try:
    from astropy.io import fits as astrofits
except ImportError:
    print("ERROR: astropy not found.")
    print("Install with: pip3 install astropy --break-system-packages")
    sys.exit(1)

hdul   = astrofits.open(fits_path)
hdu    = hdul[0]
header = hdu.header

raw = hdu.data
if   raw.ndim == 3: data = raw[0, 0, :].astype(np.float64)
elif raw.ndim == 2: data = raw[0, :].astype(np.float64)
elif raw.ndim == 1: data = raw.astype(np.float64)
else:
    print(f"ERROR: unexpected data shape {raw.shape}")
    sys.exit(1)

nchans = header['NAXIS1']
crval  = header['CRVAL1']   # Hz
cdelt  = header['CDELT1']   # Hz/channel
crpix  = header['CRPIX1']
freqs  = crval + (np.arange(nchans) - (crpix - 1)) * cdelt  # Hz

# Metadata
telescope = header.get('TELESCOP', 'unknown')
date_obs  = header.get('DATE-OBS', 'unknown')
obstime   = header.get('OBSTIME',  0)
glon      = header.get('CRVAL2',   float('nan'))
glat      = header.get('CRVAL3',   float('nan'))
bunit     = header.get('BUNIT',    '?')

restfreq = 1420405750.0       # HI rest frequency Hz
c        = 2.99792458e8       # m/s
vlsr     = (restfreq - freqs) / restfreq * c / 1000.0  # km/s

print(f"Telescope : {telescope}")
print(f"Date      : {date_obs}")
print(f"Obstime   : {obstime} s")
print(f"GLon/GLat : {glon:.1f} / {glat:.1f} deg")
print(f"Channels  : {nchans}")
print(f"Resolution: {cdelt:.2f} Hz/channel")
print(f"Freq range: {freqs[0]/1e6:.4f} - {freqs[-1]/1e6:.4f} MHz")
print(f"VLSR range: {vlsr[0]:.1f} - {vlsr[-1]:.1f} km/s")

# ── bandpass removal: running median ──────────────────────────────────────────
BP_WINDOW = 401   # channels (~49 kHz) — follows bandpass, ignores narrow features
baseline  = median_filter(data, size=BP_WINDOW, mode='nearest')
residual  = data - baseline

# ── local noise: rolling MAD ──────────────────────────────────────────────────
NOISE_WINDOW = 501  # channels (~61 kHz)
local_noise  = generic_filter(
    residual,
    lambda x: np.median(np.abs(x - np.median(x))) * 1.4826,
    size=NOISE_WINDOW, mode='nearest'
)
# Avoid division by zero at edges
local_noise  = np.where(local_noise < 1e-10, np.nan, local_noise)
snr          = residual / local_noise

global_noise = np.nanmedian(local_noise)
print(f"Global noise (median of local): {global_noise*1000:.4f} mK")

# ── narrowband peak finder ─────────────────────────────────────────────────────
EDGE      = 150   # ignore this many channels at each edge
MAX_WIDTH = 5     # maximum width in channels to be called "narrowband"

def find_narrowband(snr_arr, threshold, max_width=MAX_WIDTH, edge=EDGE):
    arr = np.where(np.isfinite(snr_arr), snr_arr, 0.0)
    arr[:edge]  = 0
    arr[-edge:] = 0
    above = np.abs(arr) > threshold
    candidates = []
    i = 0
    while i < len(arr):
        if above[i]:
            j = i
            while j < len(arr) and above[j]:
                j += 1
            width = j - i
            if width <= max_width:
                pk = i + np.argmax(np.abs(arr[i:j]))
                candidates.append((pk, width, arr[pk]))
            i = j
        else:
            i += 1
    return candidates

print()
print("Narrowband candidates (width <= 5 ch):")
print(f"{'Thresh':>8} {'Ch':>6} {'W':>3} {'Freq MHz':>13} {'VLSR km/s':>10} "
      f"{'SNR':>8} {'Noise mK':>9}")
print('-' * 65)

all_candidates = {}
for threshold in [10, 8, 6, 5]:
    cands = find_narrowband(snr, threshold)
    label = f"{threshold}sigma"
    new = [c for c in cands if c[0] not in all_candidates]
    all_candidates.update({c[0]: c for c in cands})
    print(f"  {label} : {len(cands)} total candidates")
    for ch, w, s in sorted(cands, key=lambda x: -abs(x[2])):
        freq = freqs[ch]
        v    = vlsr[ch]
        n_ch = local_noise[ch] * 1000
        print(f"  {threshold:>6}sigma {ch:>6} {w:>3} {freq/1e6:>13.6f} "
              f"{v:>10.2f} {s:>8.1f} {n_ch:>9.3f}")

# Best candidate
best = None
if all_candidates:
    best_ch, best_w, best_snr = max(all_candidates.values(), key=lambda x: abs(x[2]))
    best = (best_ch, best_w, best_snr)
    print(f"\nBest: ch {best_ch}  {freqs[best_ch]/1e6:.6f} MHz  "
          f"VLSR {vlsr[best_ch]:+.2f} km/s  SNR {best_snr:.1f}sigma  "
          f"width {best_w} ch (~{best_w*abs(cdelt):.0f} Hz)")
else:
    print("\nNo candidates above 5sigma with width <= 5 channels.")

# ── common title ──────────────────────────────────────────────────────────────
title_obs = (f"SALSA {telescope} — GLon={glon:.0f}°, GLat={glat:.0f}° — "
             f"Raw mode — {int(obstime)}s\n"
             f"{date_obs} UTC | {nchans} ch | {abs(cdelt):.0f} Hz/ch | "
             f"noise {global_noise*1000:.2f} mK (local MAD)")

# ── PLOT 1: full band ─────────────────────────────────────────────────────────
fig, axes = plt.subplots(3, 1, figsize=(14, 12))
fig.suptitle(title_obs, fontsize=12, fontweight='bold')

ax = axes[0]
ax.plot(freqs/1e6, data,     color='steelblue', lw=0.5, alpha=0.8, label='Raw spectrum')
ax.plot(freqs/1e6, baseline, color='red',       lw=1.5, ls='--',
        label=f'Bandpass (running median {BP_WINDOW} ch)')
ax.axvline(restfreq/1e6, color='orange', lw=1, ls=':', alpha=0.7, label='HI rest freq')
if best:
    ax.axvline(freqs[best[0]]/1e6, color='magenta', lw=1.2, alpha=0.7,
               label=f'Candidate {best[2]:.1f}σ')
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel(f'Amplitude ({bunit}, nominal)')
ax.set_title('Full band — raw spectrum + bandpass model')
ax.legend(fontsize=9)
ax.set_xlim(freqs[0]/1e6, freqs[-1]/1e6)
ax.grid(True, alpha=0.3)

ax = axes[1]
ax.plot(freqs/1e6, residual*1000, color='darkgreen', lw=0.5, alpha=0.8, label='Residual')
ax.plot(freqs/1e6,  local_noise*1000*5,  color='orange', lw=0.8, ls='--', alpha=0.7, label='±5σ (local)')
ax.plot(freqs/1e6, -local_noise*1000*5,  color='orange', lw=0.8, ls='--', alpha=0.7)
ax.plot(freqs/1e6,  local_noise*1000*10, color='red',    lw=0.8, ls='--', alpha=0.7, label='±10σ (local)')
ax.plot(freqs/1e6, -local_noise*1000*10, color='red',    lw=0.8, ls='--', alpha=0.7)
ax.axhline(0, color='black', lw=0.8)
if best:
    ax.axvline(freqs[best[0]]/1e6, color='magenta', lw=1.5, alpha=0.8,
               label=f'Candidate {best[2]:.1f}σ')
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Residual (mK)')
ax.set_title('Residual after bandpass subtraction — local ±5σ/±10σ envelope')
ax.legend(fontsize=9)
ax.set_xlim(freqs[0]/1e6, freqs[-1]/1e6)
ax.grid(True, alpha=0.3)

ax = axes[2]
snr_plot = np.where(np.isfinite(snr), snr, 0)
ax.plot(freqs/1e6, snr_plot, color='navy', lw=0.5, alpha=0.8)
ax.axhline(0,  color='black',  lw=0.8)
ax.axhline( 5, color='orange', lw=1, ls='--', alpha=0.7, label='±5σ')
ax.axhline(-5, color='orange', lw=1, ls='--', alpha=0.7)
ax.axhline( 10, color='red',   lw=1, ls='--', alpha=0.7, label='±10σ')
ax.axhline(-10, color='red',   lw=1, ls='--', alpha=0.7)
if best:
    ax.scatter([freqs[best[0]]/1e6], [best[2]], color='magenta', s=60,
               zorder=5, label=f'ch{best[0]}: {best[2]:.1f}σ')
ax.set_xlabel('Frequency (MHz)')
ax.set_ylabel('Local SNR (σ)')
ax.set_title('Local signal-to-noise ratio (residual / local MAD noise)')
ax.legend(fontsize=9)
ax.set_xlim(freqs[0]/1e6, freqs[-1]/1e6)
ax.set_ylim(-15, max(30, best[2]+5) if best else 15)
ax.grid(True, alpha=0.3)

plt.tight_layout()
out1 = os.path.join(plots_dir, 'SETI_analysis_full.png')
plt.savefig(out1, dpi=150, bbox_inches='tight')
plt.close()
print(f"\nSaved: {out1}")

# ── PLOT 2: zoom on best candidate ────────────────────────────────────────────
if best:
    best_ch, best_w, best_snr = best
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(
        f"Narrowband candidate — ch {best_ch} — "
        f"{freqs[best_ch]/1e6:.6f} MHz — "
        f"VLSR {vlsr[best_ch]:+.2f} km/s — SNR {best_snr:.1f}σ",
        fontsize=12, fontweight='bold')

    zoom1, zoom2 = 50, 10
    lo1 = max(0,      best_ch - zoom1)
    hi1 = min(nchans, best_ch + zoom1)
    lo2 = max(0,      best_ch - zoom2)
    hi2 = min(nchans, best_ch + zoom2)

    ax = axes[0]
    ax.plot(freqs[lo1:hi1]/1e6, residual[lo1:hi1]*1000, color='darkgreen', lw=1.2)
    ax.fill_between(freqs[lo1:hi1]/1e6,
                    -local_noise[lo1:hi1]*1000*5,
                     local_noise[lo1:hi1]*1000*5,
                    color='orange', alpha=0.15, label='±5σ envelope')
    ax.fill_between(freqs[lo1:hi1]/1e6,
                    -local_noise[lo1:hi1]*1000*10,
                     local_noise[lo1:hi1]*1000*10,
                    color='red', alpha=0.10, label='±10σ envelope')
    ax.axhline(0, color='black', lw=0.8)
    ax.axvline(freqs[best_ch]/1e6, color='magenta', lw=1.5, ls=':', label=f'ch {best_ch}')
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Residual (mK)')
    ax.set_title(f'Residual ±{zoom1} channels around candidate')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    ax = axes[1]
    colors = ['magenta' if abs(s) > 10 else 'steelblue'
              for s in snr_plot[lo2:hi2]]
    ax.bar(freqs[lo2:hi2]/1e6, snr_plot[lo2:hi2],
           width=abs(cdelt)/1e6*0.8, color=colors,
           alpha=0.85, edgecolor='black', lw=0.5)
    ax.axhline( 5,  color='orange', lw=1.5, ls='--', label='5σ')
    ax.axhline(-5,  color='orange', lw=1.5, ls='--')
    ax.axhline( 10, color='red',    lw=1.5, ls='--', label='10σ')
    ax.axhline(-10, color='red',    lw=1.5, ls='--')
    ax.set_xlabel('Frequency (MHz)')
    ax.set_ylabel('Local SNR (σ)')
    ax.set_title(f'Ultra-zoom ±{zoom2} channels — {abs(cdelt):.0f} Hz/ch')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    offset = 3 * abs(cdelt) / 1e6
    ax.annotate(
        f'{best_snr:.1f}σ\n{freqs[best_ch]/1e6:.6f} MHz\n'
        f'{vlsr[best_ch]:+.2f} km/s LSR\nwidth: {best_w} ch (~{best_w*abs(cdelt):.0f} Hz)',
        xy=(freqs[best_ch]/1e6, best_snr),
        xytext=(freqs[best_ch]/1e6 + offset, best_snr * 0.75),
        fontsize=8, color='darkred',
        arrowprops=dict(arrowstyle='->', color='darkred'))

    plt.tight_layout()
    out2 = os.path.join(plots_dir, 'SETI_candidate_zoom.png')
    plt.savefig(out2, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"Saved: {out2}")

# ── PLOT 3: SNR histogram ─────────────────────────────────────────────────────
from scipy import stats as sp_stats

fig, ax = plt.subplots(figsize=(10, 5))
snr_inner = snr_plot[EDGE:-EDGE]
snr_max   = max(25, best[2]+3) if best else 15
bins = np.linspace(-15, snr_max, 250)
ax.hist(snr_inner, bins=bins, color='steelblue', alpha=0.7,
        edgecolor='none', label='All channels')
x_g   = np.linspace(-15, 15, 500)
gauss = sp_stats.norm.pdf(x_g, 0, 1) * len(snr_inner) * (bins[1]-bins[0])
ax.plot(x_g, gauss, 'r-', lw=2, label='Expected Gaussian noise')
ax.axvline( 5, color='orange', lw=1.5, ls='--', label='5σ')
ax.axvline(10, color='red',    lw=1.5, ls='--', label='10σ')
if best:
    ax.axvline(best[2], color='magenta', lw=2,
               label=f'Candidate: {best[2]:.1f}σ')
ax.set_xlabel('Local SNR (σ)')
ax.set_ylabel('Number of channels')
ax.set_title(
    f'SNR distribution — SALSA Raw {int(obstime)}s — '
    f'GLon={glon:.0f}°, GLat={glat:.0f}°\n'
    f'Comparison with expected Gaussian noise')
ax.legend(fontsize=9)
ax.set_yscale('log')
ax.grid(True, alpha=0.3)
ax.set_xlim(-15, snr_max)

plt.tight_layout()
out3 = os.path.join(plots_dir, 'SETI_snr_histogram.png')
plt.savefig(out3, dpi=150, bbox_inches='tight')
plt.close()
print(f"Saved: {out3}")

print("\nDone.")
