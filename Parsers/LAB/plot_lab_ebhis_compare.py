# plot_lab_ebhis_compare.py
# Side-by-side LV diagram comparison: EBHIS vs LAB at b=0
# FRO/HALO Project — Francesco Di Giovanni / Claude AI (Anthropic)

import os, glob, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# ── Auto-detect available survey grids ─────────────────────────────────────────
_base = os.path.expanduser('~/FRO/Parsers/LAB/Surveys/')
_grids = sorted([d for d in os.listdir(_base) if d.startswith('LAB_Survey_G')])
if not _grids:
    raise RuntimeError('No LAB survey directory found in ' + _base)
if len(_grids) == 1:
    GRID_TAG = _grids[0].replace('LAB_Survey_', '')
else:
    print('Available grids:')
    for i, g in enumerate(_grids):
        print(f'  {i+1}. {g}')
    choice = int(input('Choose [1-N]: ').strip()) - 1
    GRID_TAG = _grids[choice].replace('LAB_Survey_', '')
SURVEY_DIR = os.path.join(_base, f'LAB_Survey_{GRID_TAG}')
LAB_DIR    = os.path.join(SURVEY_DIR, f'LAB_Profiles_{GRID_TAG}') + '/'
OUTPUT     = os.path.join(SURVEY_DIR, f'LAB_Plots_{GRID_TAG}', 'LAB_EBHIS_compare_lv.png')
# Auto-detect EBHIS grid
_ebase = os.path.expanduser('~/FRO/Parsers/Effelsberg/Surveys/')
_egrids = sorted([d for d in os.listdir(_ebase) if d.startswith('EBHIS_Survey_G')])
if not _egrids:
    raise RuntimeError('No EBHIS survey directory found')
if len(_egrids) == 1:
    _etag = _egrids[0].replace('EBHIS_Survey_', '')
else:
    print('Available EBHIS grids:')
    for i, g in enumerate(_egrids): print(f'  {i+1}. {g}')
    _etag = _egrids[int(input('Choose [1-N]: ').strip())-1].replace('EBHIS_Survey_', '')
EBHIS_DIR = os.path.join(_ebase, f'EBHIS_Survey_{_etag}', f'EBHIS_Profiles_{_etag}') + '/' 
VEL_MIN    = -120.0
VEL_MAX    =   80.0
# Extract GLON/GLAT bounds from GRID_TAG
import re as _re
_m = _re.match(r"G([\d.]+)-([\d.]+)_B([+-]?[\d.]+)([+-][\d.]+)", GRID_TAG)
GLON_MIN=float(_m.group(1)); GLON_MAX=float(_m.group(2))
GLAT_MIN=float(_m.group(3)); GLAT_MAX=float(_m.group(4))
# ───────────────────────────────────────────────────────────────────────────────

def read_ebhis_profile(filepath):
    vels, tbs = [], []
    in_block = False
    with open(filepath) as f:
        for line in f:
            if '%%EBHIS' in line:
                in_block = True
                continue
            if in_block:
                if line.startswith('%%'):
                    break
                line = line.strip()
                if line.startswith('%') or not line:
                    continue
                p = line.split()
                if len(p) >= 2:
                    vels.append(float(p[0]))
                    tbs.append(float(p[1]))
    return np.array(vels), np.array(tbs)


def read_lab_profile(filepath):
    vels, tbs = [], []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith('%') or not line:
                continue
            p = line.split()
            if len(p) == 2:
                try:
                    vels.append(float(p[0]))
                    tbs.append(float(p[1]))
                except:
                    pass
    return np.array(vels), np.array(tbs)


def extract_coords_ebhis(fname):
    m = re.match(r'EBHIS_G(\d+\.\d+)([+-]\d+\.\d+)\.txt', fname)
    if m:
        return float(m.group(1)), float(m.group(2))
    return None, None


def extract_coords_lab(fname):
    core = fname.replace('LAB_G', '').replace('.txt', '')
    try:
        if '+' in core[1:]:
            idx = core.index('+', 1)
            return float(core[:idx]), float(core[idx:])
        else:
            parts = core.split('-')
            return float(parts[0]), -float(parts[1])
    except:
        return None, None


def build_lv_grid(profiles):
    glon_bins = np.arange(GLON_MIN, GLON_MAX + 0.1, 0.5)
    vel_bins  = np.linspace(VEL_MIN, VEL_MAX, 500)
    grid = np.full((len(vel_bins), len(glon_bins)), np.nan)
    for glon, vel, tb in profiles:
        gi = np.argmin(np.abs(glon_bins - glon))
        grid[:, gi] = np.interp(vel_bins, vel, tb)
    for vi in range(grid.shape[0]):
        row = grid[vi, :]
        valid = np.where(~np.isnan(row))[0]
        if len(valid) < 2:
            continue
        fi = interp1d(valid, row[valid], kind='linear',
                      bounds_error=False, fill_value=np.nan)
        grid[vi, :] = fi(np.arange(len(glon_bins)))
    return np.clip(grid, 0, None), glon_bins, vel_bins


# Load EBHIS b=0
ebhis_profiles = []
for fp in sorted(glob.glob(os.path.join(EBHIS_DIR, 'EBHIS_G*.txt'))):
    glon, glat = extract_coords_ebhis(os.path.basename(fp))
    if glon is None or glat != 0.0:
        continue
    v, t = read_ebhis_profile(fp)
    if len(v):
        ebhis_profiles.append((glon, v, t))

# Load LAB b=0
lab_profiles = []
for fp in sorted(glob.glob(os.path.join(LAB_DIR, 'LAB_G*.txt'))):
    glon, glat = extract_coords_lab(os.path.basename(fp))
    if glon is None or glat != 0.0:
        continue
    v, t = read_lab_profile(fp)
    if len(v):
        lab_profiles.append((glon, v, t))

print(f'EBHIS profiles at b=0: {len(ebhis_profiles)}')
print(f'LAB   profiles at b=0: {len(lab_profiles)}')

ebhis_grid, glon_bins, vel_bins = build_lv_grid(ebhis_profiles)
lab_grid, _, _ = build_lv_grid(lab_profiles)

vmax = max(np.nanmax(ebhis_grid), np.nanmax(lab_grid))


while True:
    print('\nLAB vs EBHIS comparison — options:')
    print('  1. Generate LV comparison diagram')
    print('  2. Exit')
    choice = input('Choose [1-2]: ').strip()
    if choice == '2':
        print('Bye.')
        break
    if choice == '1':
        os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
        if os.path.exists(OUTPUT):
            ow = input(f'  {os.path.basename(OUTPUT)} already exists. Overwrite? [y/N]: ').strip().lower()
            if ow != 'y':
                print('  Skipped.')
                continue
        fig, axes = plt.subplots(2, 1, figsize=(14, 10), sharex=True, sharey=True)
        for ax, grid, title in zip(axes,
            [ebhis_grid, lab_grid],
            ["EBHIS (Effelsberg 100m, beam 10.8') — b=0°",
             "LAB (Dwingeloo 25m / Villa Elisa 30m, beam 36') — b=0°"]):
            im = ax.imshow(grid, origin='lower', aspect='auto',
                           extent=[GLON_MIN, GLON_MAX, VEL_MIN, VEL_MAX],
                           cmap='gnuplot', vmin=0, vmax=vmax)
            plt.colorbar(im, ax=ax, label='T_B (K)')
            ax.set_ylabel('LSR velocity (km/s)')
            ax.set_title(title)
            ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
        axes[1].set_xlabel('Galactic longitude l (deg)')
        fig.suptitle(f'HI survey comparison: EBHIS vs LAB — GLon {GLON_MIN:.0f}°–{GLON_MAX:.0f}°, b=0°', fontsize=13)
        plt.tight_layout()
        plt.savefig(OUTPUT, dpi=150)
        plt.show()
        print('Saved:', OUTPUT)
