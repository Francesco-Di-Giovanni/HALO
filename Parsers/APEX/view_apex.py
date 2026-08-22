#!/usr/bin/env python3
"""
view_apex.py — APEX MBFITS 12CO(3-2) visualizer
FRO/HALO Project — Francesco Di Giovanni (Claude AI assisted)

Produces 5 plots from APEX bb1+bb2 .fro files:
  1. Combined bb1+bb2 mean spectrum (full band ~4 GHz)
  2. bb2 mean spectrum zoomed on 12CO(3-2) line
  3. GLat-velocity (b-v) diagram — bb2
  4. GLon-velocity (l-v) diagram — bb2
  5. Strip map: GLon vs GLat colored by peak temperature

Usage: python3 view_apex.py
"""

import os
import sys
import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import subprocess

# ── Constants ──────────────────────────────────────────────────────────────────
REST_FREQ_HZ = 345.7959899e9   # 12CO(3-2) rest frequency (Hz)
C_MPS        = 299792458.0     # speed of light (m/s)

# ── Paths ───────────────────────────────────────────────────────────────────────
SCRIPT_DIR  = os.path.dirname(os.path.abspath(__file__))
SURVEY_DIR  = os.path.join(SCRIPT_DIR, 'Surveys', 'MBFITS_Obs_12CO_345GHz')
FRO_DIR     = os.path.join(SURVEY_DIR, 'MBFITS_FRO_12CO_345GHz')
PLOTS_DIR   = os.path.join(SURVEY_DIR, 'MBFITS_Plots_12CO_345GHz')

os.makedirs(PLOTS_DIR, exist_ok=True)

def safe_path(path):
    """Raise if output file already exists to prevent silent overwrite."""
    if os.path.exists(path):
        raise FileExistsError(
            f"Output file already exists: {path}\n"
            f"Remove or rename it before re-running."
        )
    return path

# ── Find .fro files ─────────────────────────────────────────────────────────────
fro_files = sorted([
    os.path.join(FRO_DIR, f)
    for f in os.listdir(FRO_DIR)
    if f.endswith('.fro')
])
assert len(fro_files) == 2, f"Expected 2 .fro files (bb1, bb2), found {len(fro_files)}"

bb1_path = [f for f in fro_files if 'bb1' in f][0]
bb2_path = [f for f in fro_files if 'bb2' in f][0]

print(f"bb1: {os.path.basename(bb1_path)}")
print(f"bb2: {os.path.basename(bb2_path)}")

# ── Load data ───────────────────────────────────────────────────────────────────
def load_fro(path):
    with h5py.File(path, 'r') as f:
        freq = f['Spectra/frequency_axis_hz'][:]
        spec = f['Spectra/raw_spectrum'][:]
        glon = f['Pointing/glon_deg'][:]
        glat = f['Pointing/glat_deg'][:]
        utc  = f['Time/utc'][:]
        facility   = f['Source'].attrs.get('facility',   'APEX')
        instrument = f['Source'].attrs.get('instrument', 'HET345-XFFTS2')
        dataset    = f['Source'].attrs.get('dataset',    '')
    return freq, spec, glon, glat, utc, facility, instrument, dataset

freq1, spec1, glon1, glat1, utc1, facility, instrument, dataset = load_fro(bb1_path)
freq2, spec2, glon2, glat2, utc2, _, _, _ = load_fro(bb2_path)

assert spec1.shape[0] == spec2.shape[0], "bb1 and bb2 must have same number of spectra"
N = spec1.shape[0]

# ── Velocity axis (LSRK) for bb2 ────────────────────────────────────────────────
vel2_kms = (REST_FREQ_HZ - freq2) / REST_FREQ_HZ * C_MPS / 1e3

# ── Common plot style ────────────────────────────────────────────────────────────
TITLE_PREFIX = f"APEX {instrument} — 12CO(3-2) — {dataset}"
CALNOTE      = "Antenna temperature (K, nominal Tsys ~300 K, not calibrated)"

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.titlesize': 11,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
})

generated = []

# ══════════════════════════════════════════════════════════════════════════════
# Plot 1 — Combined bb1+bb2 mean spectrum (full ~4 GHz band)
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 1: combined mean spectrum bb1+bb2...")

mean1 = np.nanmean(spec1, axis=0)
mean2 = np.nanmean(spec2, axis=0)

fig, ax = plt.subplots(figsize=(12, 4))
ax.plot(freq1 / 1e9, mean1, lw=0.6, color='steelblue',  label='bb1')
ax.plot(freq2 / 1e9, mean2, lw=0.6, color='darkorange', label='bb2')
ax.axvline(REST_FREQ_HZ / 1e9, color='red', lw=1.0, ls='--', label=f'12CO(3-2) {REST_FREQ_HZ/1e9:.4f} GHz')
ax.set_xlabel('Frequency (GHz)')
ax.set_ylabel('Antenna temperature (K)')
ax.set_title(f"{TITLE_PREFIX}\nMean spectrum — full band (bb1 + bb2 combined)")
ax.legend(fontsize=8)
fig.text(0.5, -0.02, CALNOTE, ha='center', fontsize=7, color='gray')
fig.tight_layout()

out1_svg = safe_path(os.path.join(PLOTS_DIR, 'APEX_12CO_mean_spectrum_full_band.svg'))
fig.savefig(out1_svg, bbox_inches='tight')
plt.close(fig)
generated.append(out1_svg)
print(f"  -> {out1_svg}")

# ══════════════════════════════════════════════════════════════════════════════
# Plot 2 — bb2 mean spectrum zoomed on 12CO(3-2) line
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 2: bb2 mean spectrum zoom on 12CO(3-2)...")

mean2_zoom = np.nanmean(spec2, axis=0)

mask_zoom = np.abs(vel2_kms) < 200.0
vel_zoom  = vel2_kms[mask_zoom]
spec_zoom = mean2_zoom[mask_zoom]

fig, ax = plt.subplots(figsize=(9, 4))
ax.plot(vel_zoom, spec_zoom, lw=0.8, color='darkorange')
ax.axvline(0.0, color='red', lw=1.0, ls='--', label='v_LSRK = 0 km/s')
ax.fill_between(vel_zoom, spec_zoom, alpha=0.15, color='darkorange')
ax.set_xlabel('LSRK velocity (km/s)')
ax.set_ylabel('Antenna temperature (K)')
ax.set_title(f"{TITLE_PREFIX}\nMean spectrum — zoom +-200 km/s around 12CO(3-2)")
ax.legend(fontsize=8)
fig.text(0.5, -0.02, CALNOTE, ha='center', fontsize=7, color='gray')
fig.tight_layout()

out2_svg = safe_path(os.path.join(PLOTS_DIR, 'APEX_12CO_mean_spectrum_zoom.svg'))
fig.savefig(out2_svg, bbox_inches='tight')
plt.close(fig)
generated.append(out2_svg)
print(f"  -> {out2_svg}")

# ══════════════════════════════════════════════════════════════════════════════
# Plot 3 — GLat-velocity (b-v) diagram — bb2
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 3: GLat-velocity diagram...")

sort_idx = np.argsort(glat2)
spec2_sorted = spec2[sort_idx, :]
glat_sorted  = glat2[sort_idx]

mask_v = np.abs(vel2_kms) < 200.0
spec2_bv = spec2_sorted[:, mask_v]
vel_bv   = vel2_kms[mask_v]

vmin = np.nanpercentile(spec2_bv, 2)
vmax = np.nanpercentile(spec2_bv, 98)

fig, ax = plt.subplots(figsize=(9, 6))
im = ax.imshow(
    spec2_bv,
    origin='lower',
    aspect='auto',
    extent=[vel_bv.min(), vel_bv.max(), glat_sorted.min(), glat_sorted.max()],
    cmap='inferno',
    vmin=vmin, vmax=vmax,
    interpolation='nearest'
)
ax.axvline(0.0, color='cyan', lw=0.8, ls='--', alpha=0.7, label='v=0 km/s')
ax.set_xlabel('LSRK velocity (km/s)')
ax.set_ylabel('Galactic latitude b (deg)')
ax.set_title(f"{TITLE_PREFIX}\nGLat-velocity diagram (bb2)")
cb = fig.colorbar(im, ax=ax, pad=0.02)
cb.set_label('Antenna temperature (K)')
ax.legend(fontsize=8)
fig.text(0.5, -0.02, CALNOTE, ha='center', fontsize=7, color='gray')
fig.tight_layout()

out3 = safe_path(os.path.join(PLOTS_DIR, 'APEX_12CO_bv_diagram.png'))
fig.savefig(out3, dpi=300, bbox_inches='tight')
plt.close(fig)
generated.append(out3)
print(f"  -> {out3}")

# ══════════════════════════════════════════════════════════════════════════════
# Plot 4 — GLon-velocity (l-v) diagram — bb2
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 4: GLon-velocity diagram...")

sort_idx_l  = np.argsort(glon2)
spec2_lv    = spec2[sort_idx_l, :][:, mask_v]
glon_sorted = glon2[sort_idx_l]

vmin_lv = np.nanpercentile(spec2_lv, 2)
vmax_lv = np.nanpercentile(spec2_lv, 98)

fig, ax = plt.subplots(figsize=(9, 6))
im = ax.imshow(
    spec2_lv,
    origin='lower',
    aspect='auto',
    extent=[vel_bv.min(), vel_bv.max(), glon_sorted.min(), glon_sorted.max()],
    cmap='inferno',
    vmin=vmin_lv, vmax=vmax_lv,
    interpolation='nearest'
)
ax.axvline(0.0, color='cyan', lw=0.8, ls='--', alpha=0.7, label='v=0 km/s')
ax.set_xlabel('LSRK velocity (km/s)')
ax.set_ylabel('Galactic longitude l (deg)')
ax.set_title(f"{TITLE_PREFIX}\nGLon-velocity diagram (bb2)")
cb = fig.colorbar(im, ax=ax, pad=0.02)
cb.set_label('Antenna temperature (K)')
ax.legend(fontsize=8)
fig.text(0.5, -0.02, CALNOTE, ha='center', fontsize=7, color='gray')
fig.tight_layout()

out4 = safe_path(os.path.join(PLOTS_DIR, 'APEX_12CO_lv_diagram.png'))
fig.savefig(out4, dpi=300, bbox_inches='tight')
plt.close(fig)
generated.append(out4)
print(f"  -> {out4}")

# ══════════════════════════════════════════════════════════════════════════════
# Plot 5 — Strip map: GLon vs GLat colored by peak temperature (bb2)
# ══════════════════════════════════════════════════════════════════════════════
print("Plot 5: strip map GLon vs GLat colored by peak temperature...")

peak_T = np.nanmax(spec2[:, mask_v], axis=1)

fig, ax = plt.subplots(figsize=(7, 6))
sc = ax.scatter(
    glon2, glat2,
    c=peak_T,
    cmap='inferno',
    s=60,
    edgecolors='none',
    vmin=np.nanpercentile(peak_T, 2),
    vmax=np.nanpercentile(peak_T, 98)
)
cb = fig.colorbar(sc, ax=ax, pad=0.02)
cb.set_label('Peak antenna temperature (K)')
ax.set_xlabel('Galactic longitude l (deg)')
ax.set_ylabel('Galactic latitude b (deg)')
ax.set_title(f"{TITLE_PREFIX}\nStrip map — peak temperature per pointing (bb2)")
ax.invert_xaxis()
fig.text(0.5, -0.02, CALNOTE, ha='center', fontsize=7, color='gray')
fig.tight_layout()

out5 = safe_path(os.path.join(PLOTS_DIR, 'APEX_12CO_strip_map_peak_T.png'))
fig.savefig(out5, dpi=300, bbox_inches='tight')
plt.close(fig)
generated.append(out5)
print(f"  -> {out5}")

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
print("\n── Plot generati ──────────────────────────────────────────────────────")
for p in generated:
    print(f"  {p}")
print("───────────────────────────────────────────────────────────────────────")

for p in generated:
    subprocess.Popen(['xdg-open', p])
