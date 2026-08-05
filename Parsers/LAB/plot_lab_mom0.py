# plot_lab_mom0.py
# HI maps from LAB survey profiles: moment-0, peak temperature, peak velocity
# FRO/HALO Project — Francesco Di Giovanni / Claude AI (Anthropic)

import os, glob
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import zoom, generic_filter

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
LAB_DIR    = os.path.join(SURVEY_DIR, f'LAB_Profiles_{GRID_TAG}')
OUTPUT_DIR = os.path.join(SURVEY_DIR, f'LAB_Plots_{GRID_TAG}')
VEL_MIN    = -120.0
VEL_MAX    =   80.0
# Extract GLON/GLAT bounds from GRID_TAG
import re as _re
_m = _re.match(r"G([\d.]+)-([\d.]+)_B([+-]?[\d.]+)([+-][\d.]+)", GRID_TAG)
GLON_MIN=float(_m.group(1)); GLON_MAX=float(_m.group(2))
GLAT_MIN=float(_m.group(3)); GLAT_MAX=float(_m.group(4))
# ───────────────────────────────────────────────────────────────────────────────

def extract_coords(fname):
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


def read_profile(filepath):
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


def fill_nan(arr):
    mask = np.isnan(arr)
    arr[mask] = generic_filter(arr, np.nanmean, size=3)[mask]
    return arr


def smooth(grid):
    grid = fill_nan(grid.copy())
    return zoom(grid, 8, order=1)


def plot_map(grid_s, glon_u, glat_u, title, label, cmap, output):
    fig, ax = plt.subplots(figsize=(14, 8))
    im = ax.imshow(grid_s, origin='lower', aspect='auto',
                   extent=[glon_u[0], glon_u[-1], glat_u[0], glat_u[-1]],
                   cmap=cmap)
    plt.colorbar(im, ax=ax, label=label)
    ax.set_xlabel('Galactic longitude l (deg)')
    ax.set_ylabel('Galactic latitude b (deg)')
    ax.set_title(title + '\nDwingeloo 25m / IAR Villa Elisa 30m — Kalberla et al. 2005\n'
                 f'GLon {GLON_MIN:.0f}°–{GLON_MAX:.0f}°, GLat {GLAT_MIN:.0f}°/+{GLAT_MAX:.0f}°, '
                 f'v = {VEL_MIN} to {VEL_MAX} km/s')
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    plt.tight_layout()
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    plt.savefig(output, dpi=150)
    print('Saved:', output)
    plt.show()


# Load all profiles
glons_u, glats_u = [], []
mom0_data, tpeak_data, vpeak_data = {}, {}, {}

for fp in sorted(glob.glob(os.path.join(LAB_DIR, 'LAB_G*.txt'))):
    glon, glat = extract_coords(os.path.basename(fp))
    if glon is None:
        continue
    v, t = read_profile(fp)
    if len(v) == 0:
        continue
    mask = (v >= VEL_MIN) & (v <= VEL_MAX)
    if mask.sum() < 2:
        continue
    dv = np.abs(np.mean(np.diff(v[mask])))
    mom0_data[(glon, glat)]  = np.sum(t[mask]) * dv
    tpeak_data[(glon, glat)] = np.max(t[mask])
    vpeak_data[(glon, glat)] = v[mask][np.argmax(t[mask])]
    glons_u.append(glon)
    glats_u.append(glat)

glons_u = sorted(set(glons_u))
glats_u = sorted(set(glats_u))


def build_grid(data_dict):
    g = np.full((len(glats_u), len(glons_u)), np.nan)
    for (glon, glat), val in data_dict.items():
        gi = glons_u.index(glon)
        bi = glats_u.index(glat)
        g[bi, gi] = val
    return g


mom0_grid  = build_grid(mom0_data)
tpeak_grid = build_grid(tpeak_data)
vpeak_grid = build_grid(vpeak_data)

# Menu — ciclic until Exit
while True:
    print('\nLAB survey — available maps:')
    print('  1. Moment-0 map (integrated intensity, K km/s)')
    print('  2. Peak temperature map (K)')
    print('  3. Peak velocity map (km/s)')
    print('  4. All three (separate figures)')
    print('  5. All three (single figure)')
    print('  6. Exit')
    choice = input('Choose [1-6]: ').strip()

    if choice == '6':
        print('Bye.')
        break

    def _out(name):
        p = os.path.join(OUTPUT_DIR, name)
        if os.path.exists(p):
            ow = input(f'  {name} already exists. Overwrite? [y/N]: ').strip().lower()
            if ow != 'y':
                print('  Skipped.')
                return None
        return p

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    if choice in ('1', '4'):
        out = _out('LAB_mom0_map.png')
        if out: plot_map(smooth(mom0_grid), glons_u, glats_u,
                 'LAB survey — HI moment-0 map',
                 'Integrated intensity (K km/s)', 'inferno', out)

    if choice in ('2', '4'):
        out = _out('LAB_tpeak_map.png')
        if out: plot_map(smooth(tpeak_grid), glons_u, glats_u,
                 'LAB survey — HI peak temperature map',
                 'Peak T_B (K)', 'hot', out)

    if choice in ('3', '4'):
        out = _out('LAB_vpeak_map.png')
        if out: plot_map(smooth(vpeak_grid), glons_u, glats_u,
                 'LAB survey — HI peak velocity map',
                 'Peak velocity (km/s)', 'RdBu_r', out)

    if choice == '5':
        out = _out('LAB_all_maps.png')
        if out:
            maps = [
                (smooth(mom0_grid),  'HI moment-0',         'Integrated intensity (K km/s)', 'inferno'),
                (smooth(tpeak_grid), 'HI peak temperature', 'Peak T_B (K)',                  'hot'),
                (smooth(vpeak_grid), 'HI peak velocity',    'Peak velocity (km/s)',           'RdBu_r'),
            ]
            fig, axes = plt.subplots(3, 1, figsize=(14, 18), sharex=True)
            for ax, (grid_s, title, label, cmap) in zip(axes, maps):
                im = ax.imshow(grid_s, origin='lower', aspect='auto',
                               extent=[glons_u[0], glons_u[-1], glats_u[0], glats_u[-1]],
                               cmap=cmap)
                plt.colorbar(im, ax=ax, label=label)
                ax.set_ylabel('Galactic latitude b (deg)')
                ax.set_title(f'LAB survey — {title}')
                ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
            axes[-1].set_xlabel('Galactic longitude l (deg)', labelpad=10)
            fig.suptitle('LAB survey — HI maps\nDwingeloo 25m / IAR Villa Elisa 30m — Kalberla et al. 2005\n'
                         f'GLon {GLON_MIN:.0f}°–{GLON_MAX:.0f}°, GLat {GLAT_MIN:.0f}°/+{GLAT_MAX:.0f}°, '
                         f'v = {VEL_MIN} to {VEL_MAX} km/s', fontsize=12)
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            plt.savefig(out, dpi=150)
            print('Saved:', out)
            plt.show()