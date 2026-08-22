#!/usr/bin/env python3
"""
view_apex_interactive.py — APEX MBFITS 12CO(3-2) interactive spectrum viewer
FRO/HALO Project — Francesco Di Giovanni (Claude AI assisted)

Two interactive plots with crosshair cursor:
  1. Combined bb1+bb2 mean spectrum (full band ~4 GHz) — frequency axis
  2. bb2 mean spectrum zoomed on 12CO(3-2) — velocity axis

Usage: python3 view_apex_interactive.py
"""

import os
import h5py
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

# ── Constants ──────────────────────────────────────────────────────────────────
REST_FREQ_HZ = 345.7959899e9
C_MPS        = 299792458.0

# ── Paths ───────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FRO_DIR    = os.path.join(SCRIPT_DIR, 'Surveys', 'MBFITS_Obs_12CO_345GHz', 'MBFITS_FRO_12CO_345GHz')

fro_files = sorted([os.path.join(FRO_DIR, f) for f in os.listdir(FRO_DIR) if f.endswith('.fro')])
bb1_path = [f for f in fro_files if 'bb1' in f][0]
bb2_path = [f for f in fro_files if 'bb2' in f][0]

def load_fro(path):
    with h5py.File(path, 'r') as f:
        freq = f['Spectra/frequency_axis_hz'][:]
        spec = f['Spectra/raw_spectrum'][:]
    return freq, spec

freq1, spec1 = load_fro(bb1_path)
freq2, spec2 = load_fro(bb2_path)
mean1 = np.nanmean(spec1, axis=0)
mean2 = np.nanmean(spec2, axis=0)
vel2_kms = (REST_FREQ_HZ - freq2) / REST_FREQ_HZ * C_MPS / 1e3
mask_zoom = np.abs(vel2_kms) < 200.0
vel_zoom  = vel2_kms[mask_zoom]
spec_zoom = mean2[mask_zoom]

CALNOTE = "Antenna temperature (K, nominal Tsys ~300 K, not calibrated)"

# ── Plot 1 — full band ─────────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(13, 5))
fig1.canvas.manager.set_window_title('APEX — Full band spectrum (bb1+bb2)')
ax1.set_xlim(343.79, 347.80)
ax1.autoscale(enable=False, axis='x')
l1a, = ax1.plot(freq1/1e9, mean1, lw=0.7, color='steelblue',  label='bb1')
l1b, = ax1.plot(freq2/1e9, mean2, lw=0.7, color='darkorange', label='bb2')
ax1.axvline(REST_FREQ_HZ/1e9, color='red', lw=1.0, ls='--', label='12CO(3-2)')
ax1.set_xlabel('Frequency (GHz)')
ax1.set_ylabel('Antenna temperature (K)')
ax1.set_title('APEX HET345-XFFTS2 — 12CO(3-2)\nMean spectrum — full band (bb1 + bb2)')
ax1.legend(fontsize=8)
fig1.text(0.5, 0.01, CALNOTE, ha='center', fontsize=7, color='gray')

vline1 = ax1.axvline(color='gray', lw=0.8, ls=':', alpha=0.7)
hline1 = ax1.axhline(color='gray', lw=0.8, ls=':', alpha=0.7)
text1  = ax1.text(0.02, 0.95, '', transform=ax1.transAxes,
                  fontsize=9, va='top', family='monospace',
                  bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))

def on_move1(event):
    if event.inaxes != ax1:
        return
    x = event.xdata
    y = event.ydata
    if x is None or y is None:
        return
    vline1.set_xdata([x])
    hline1.set_ydata([y])
    # find nearest bb1 or bb2 value
    idx1 = np.argmin(np.abs(freq1/1e9 - x))
    idx2 = np.argmin(np.abs(freq2/1e9 - x))
    v1 = mean1[idx1]
    v2 = mean2[idx2]
    text1.set_text(f'f = {x:.5f} GHz\nbb1 = {v1:.3e} K\nbb2 = {v2:.3e} K')
    fig1.canvas.draw_idle()

fig1.canvas.mpl_connect('motion_notify_event', on_move1)
fig1.tight_layout(rect=[0, 0.04, 1, 1])

# ── Plot 2 — zoom velocity ─────────────────────────────────────────────────────
fig2, ax2 = plt.subplots(figsize=(10, 5))
fig2.canvas.manager.set_window_title('APEX — Zoom 12CO(3-2) ±200 km/s')
ax2.plot(vel_zoom, spec_zoom, lw=0.9, color='darkorange')
ax2.axvline(0.0, color='red', lw=1.0, ls='--', label='v_LSRK = 0 km/s')
ax2.fill_between(vel_zoom, spec_zoom, alpha=0.15, color='darkorange')
ax2.set_xlabel('LSRK velocity (km/s)')
ax2.set_ylabel('Antenna temperature (K)')
ax2.set_title('APEX HET345-XFFTS2 — 12CO(3-2)\nMean spectrum — zoom ±200 km/s')
ax2.legend(fontsize=8)
fig2.text(0.5, 0.01, CALNOTE, ha='center', fontsize=7, color='gray')

vline2 = ax2.axvline(color='gray', lw=0.8, ls=':', alpha=0.7)
hline2 = ax2.axhline(color='gray', lw=0.8, ls=':', alpha=0.7)
text2  = ax2.text(0.02, 0.95, '', transform=ax2.transAxes,
                  fontsize=9, va='top', family='monospace',
                  bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.8))

def on_move2(event):
    if event.inaxes != ax2:
        return
    x = event.xdata
    y = event.ydata
    if x is None or y is None:
        return
    vline2.set_xdata([x])
    hline2.set_ydata([y])
    idx = np.argmin(np.abs(vel_zoom - x))
    T   = spec_zoom[idx]
    # corrispondente frequenza
    f_ghz = REST_FREQ_HZ * (1 - vel_zoom[idx]*1e3/C_MPS) / 1e9
    text2.set_text(f'v = {x:.2f} km/s\nf = {f_ghz:.5f} GHz\nT = {T:.3e} K')
    fig2.canvas.draw_idle()

fig2.canvas.mpl_connect('motion_notify_event', on_move2)
fig2.tight_layout(rect=[0, 0.04, 1, 1])

plt.show()
