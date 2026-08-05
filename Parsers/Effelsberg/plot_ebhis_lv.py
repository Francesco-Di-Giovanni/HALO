# plot_ebhis_lv.py
# LV diagram from EBHIS/HI4PI downloaded profiles.
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy

import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt

# ── Auto-detect available survey grids ─────────────────────────────────────────
_base = os.path.expanduser('~/FRO/Parsers/Effelsberg/Surveys/')
_grids = sorted([d for d in os.listdir(_base) if d.startswith('EBHIS_Survey_G')])
if not _grids:
    raise RuntimeError('No EBHIS survey directory found in ' + _base)
if len(_grids) == 1:
    GRID_TAG = _grids[0].replace('EBHIS_Survey_', '')
else:
    print('Available grids:')
    for i, g in enumerate(_grids):
        print(f'  {i+1}. {g}')
    choice = int(input('Choose [1-N]: ').strip()) - 1
    GRID_TAG = _grids[choice].replace('EBHIS_Survey_', '')
SURVEY_DIR = os.path.join(_base, f'EBHIS_Survey_{GRID_TAG}')
EBHIS_DIR  = os.path.join(SURVEY_DIR, f'EBHIS_Profiles_{GRID_TAG}') + '/'
OUTPUT     = os.path.join(SURVEY_DIR, f'EBHIS_Plots_{GRID_TAG}', 'EBHIS_lv_diagram.png')
VEL_MIN    = -120.0
VEL_MAX    =   80.0
# Extract GLON/GLAT bounds from GRID_TAG
import re as _re
_m = _re.match(r"G([\d.]+)-([\d.]+)_B([+-]?[\d.]+)([+-][\d.]+)", GRID_TAG)
GLON_MIN=float(_m.group(1)); GLON_MAX=float(_m.group(2))
GLAT_MIN=float(_m.group(3)); GLAT_MAX=float(_m.group(4))
# ───────────────────────────────────────────────────────────────────────────────

INTERPOLATE = True


def parse_ebhis_section(filepath):
    """Extract EBHIS velocity and T_B data from a downloaded profile file."""
    velocities = []
    temperatures = []
    in_ebhis = False

    with open(filepath, 'r') as f:
        for line in f:
            if '%%EBHIS' in line:
                in_ebhis = True
                continue
            if in_ebhis:
                if line.startswith('%%'):
                    break
                line = line.strip()
                if line.startswith('%') or not line:
                    continue
                parts = line.split()
                if len(parts) >= 2:
                    velocities.append(float(parts[0]))
                    temperatures.append(float(parts[1]))

    return np.array(velocities), np.array(temperatures)


def extract_coords(filename):
    """Extract GLon and GLat from filename like EBHIS_G080.00+05.00.txt"""
    match = re.match(r'EBHIS_G(\d+\.\d+)([+-]\d+\.\d+)\.txt', filename)
    if match:
        return float(match.group(1)), float(match.group(2))
    return None, None


def main():
    files = sorted(glob.glob(os.path.join(EBHIS_DIR, 'EBHIS_G*.txt')))
    print('Files found:', len(files))

    # Read all profiles
    profiles = []
    for filepath in files:
        glon, glat = extract_coords(os.path.basename(filepath))
        if glon is None:
            continue
        if glat != 0.0:
            continue   # only b=0 for the lv diagram
        vel, tb = parse_ebhis_section(filepath)
        if len(vel) == 0:
            print('  WARNING: no EBHIS data in', os.path.basename(filepath))
            continue
        profiles.append((glon, vel, tb))

    print('Profiles at b=0:', len(profiles))

    # Build LV grid
    glon_bins = np.arange(GLON_MIN, GLON_MAX + 0.1, 0.5)
    vel_bins  = np.linspace(VEL_MIN, VEL_MAX, 500)
    lv_grid   = np.full((len(vel_bins), len(glon_bins)), np.nan)

    for glon, vel, tb in profiles:
        gi = np.argmin(np.abs(glon_bins - glon))
        tb_interp = np.interp(vel_bins, vel, tb)
        lv_grid[:, gi] = tb_interp

    # Interpolate
    if INTERPOLATE:
        from scipy.interpolate import interp1d
        for vi in range(lv_grid.shape[0]):
            row = lv_grid[vi, :]
            valid = np.where(~np.isnan(row))[0]
            if len(valid) < 2:
                continue
            f_interp = interp1d(valid, row[valid], kind='linear',
                                bounds_error=False, fill_value=np.nan)
            lv_grid[vi, :] = f_interp(np.arange(len(glon_bins)))

    # Clip negatives
    lv_grid = np.clip(lv_grid, 0, None)

    # Plot
    fig, ax = plt.subplots(figsize=(14, 7))
    im = ax.imshow(lv_grid, origin='lower', aspect='auto',
                   extent=[GLON_MIN, GLON_MAX, VEL_MIN, VEL_MAX],
                   cmap='gnuplot')
    plt.colorbar(im, ax=ax, label='T_B (K)')
    ax.set_xlabel('Galactic longitude l (deg)')
    ax.set_ylabel('LSR velocity (km/s)')
    ax.set_title('EBHIS (Effelsberg 100m) — LV diagram at b=0\nHI4PI / Argelander-Institut Bonn')
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    plt.tight_layout()
    plt.savefig(OUTPUT, dpi=150)
    plt.show()
    print('Saved:', OUTPUT)


if __name__ == '__main__':
    main()