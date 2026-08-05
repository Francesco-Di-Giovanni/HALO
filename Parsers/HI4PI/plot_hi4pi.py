# plot_hi4pi.py
# Visualisation suite for HI4PI position-position-velocity cubes and for the
# .fro spectra extracted from them by hi4pi_to_fro.py.
#
# Nothing is hardcoded: cube path, latitude of the l-v cut and velocity
# integration range are all supplied at runtime.
#
# All figures are written to Surveys/HI4PI_ExtractedPoints/Plots/
#
# Author: Francesco Di Giovanni / HALO project
# Date: 2026-07-27

import os
import sys
import glob
import subprocess
import numpy as np
import h5py
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from astropy.io import fits
from astropy.wcs import WCS

C_MPS = 299792458.0
NHI_FACTOR = 1.823e18  # cm^-2 / (K km/s), optically thin HI

SURVEYS_ROOT = os.path.expanduser('~/FRO/Parsers/HI4PI/Surveys')

# Set at runtime from the chosen cube file: no fixed working directory.
FRO_DIR = ''
PLOT_DIR = ''

def safe_path(directory, fname):
    """Return a path that does not overwrite an existing file.
    If directory/fname exists, appends _1, _2, ... before the extension."""
    base, ext = os.path.splitext(fname)
    candidate = os.path.join(directory, fname)
    counter = 1
    while os.path.exists(candidate):
        candidate = os.path.join(directory, f"{base}_{counter}{ext}")
        counter += 1
    return candidate



def get_fits_path():
    if len(sys.argv) > 1:
        return os.path.expanduser(sys.argv[1])
    return os.path.expanduser(input('Path to HI4PI FITS cube: ').strip())


def load_cube(filepath):
    with fits.open(filepath) as hdul:
        data = hdul[0].data.astype(np.float32)
        header = hdul[0].header.copy()
    return data, header


def velocity_axis_kms(header):
    """LSRK radio velocity axis in km/s, from the header only."""
    n = header['NAXIS3']
    pix = np.arange(1, n + 1)
    v_mps = header['CRVAL3'] + (pix - header['CRPIX3']) * header['CDELT3']
    return v_mps / 1000.0


def world_axes(header, wcs2d):
    """Return the GLon and GLat values along the two spatial axes.
    Handles the 360/0 degree wrap: if the longitude axis is not monotonic
    (tile crosses l=0), values that wrapped below 180 are shifted by +360
    to restore continuity."""
    nx = header['NAXIS1']
    ny = header['NAXIS2']
    ix = np.arange(nx)
    iy = np.arange(ny)
    glon, _ = wcs2d.wcs_pix2world(ix, np.zeros(nx), 0)
    _, glat = wcs2d.wcs_pix2world(np.zeros(ny), iy, 0)
    glon = np.asarray(glon, dtype=float)
    glat = np.asarray(glat, dtype=float)
    # Fix 360/0 wrap: if longitude is not monotonic, unwrap it
    if len(glon) > 1:
        diffs = np.diff(glon)
        if not (np.all(diffs > 0) or np.all(diffs < 0)):
            glon = np.unwrap(glon, period=360.0)
    return glon, glat


def ask_velocity_range(v_kms):
    lo, hi = float(v_kms.min()), float(v_kms.max())
    txt = input(f'Velocity range km/s [blank = full {lo:.0f}..{hi:.0f}]: ').strip()
    if not txt:
        return lo, hi
    try:
        a, b = [float(x) for x in txt.replace(',', ' ').split()]
        return min(a, b), max(a, b)
    except ValueError:
        print('  Unrecognised range, using full range.')
        return lo, hi


# ---------------------------------------------------------------- plots ----

def plot_extracted_spectra(tag):
    """One figure per .fro previously produced by hi4pi_to_fro.py."""
    files = sorted(glob.glob(os.path.join(FRO_DIR, '*.fro')))
    if not files:
        print(f'  No .fro files found in {FRO_DIR}')
        return []
    written = []
    for path in files:
        with h5py.File(path, 'r') as f:
            freq = f['Spectra/frequency_axis_hz'][:]
            spec = f['Spectra/raw_spectrum'][0, :]
            glon = f['Pointing/glon_deg'][0]
            glat = f['Pointing/glat_deg'][0]
            restfreq = f['Observation'].attrs.get('center_frequency_hz', np.nan)
        if not np.isfinite(restfreq):
            restfreq = 1420405751.77
        v_kms = C_MPS * (1.0 - freq / restfreq) / 1000.0
        order = np.argsort(v_kms)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(v_kms[order], spec[order], lw=1.0, color='#1f4e9c')
        ax.set_xlabel('LSRK velocity (km/s)')
        ax.set_ylabel('Brightness temperature $T_B$ (K)')
        ax.set_title(f'HI4PI spectrum  GLon={glon:.2f}  GLat={glat:.2f}')
        ax.grid(alpha=0.3)
        fig.tight_layout()
        out = safe_path(PLOT_DIR, f'HI4PI_spectrum_G{glon:06.2f}{glat:+06.2f}.png')
        fig.savefig(out, dpi=130)
        plt.close(fig)
        written.append(out)
    return written


def plot_lv(data, header, wcs2d, v_kms, tag):
    glon_ax, glat_ax = world_axes(header, wcs2d)
    print(f'  GLat available: {glat_ax.min():.2f} .. {glat_ax.max():.2f}')
    txt = input('  GLat of the l-v cut (deg): ').strip()
    try:
        glat_req = float(txt)
    except ValueError:
        print('  Invalid latitude, skipping.')
        return []
    iy = int(np.argmin(np.abs(glat_ax - glat_req)))
    glat_actual = glat_ax[iy]

    v_lo, v_hi = ask_velocity_range(v_kms)
    sel = (v_kms >= v_lo) & (v_kms <= v_hi)
    v_cut = v_kms[sel]

    slab = data[sel, iy, :]  # (velocity, glon)

    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        slab, origin='lower', aspect='auto', cmap='inferno',
        extent=[glon_ax[0], glon_ax[-1], v_cut[0], v_cut[-1]])
    ax.set_xlabel('Galactic longitude $l$ (deg)')
    ax.set_ylabel('LSRK velocity (km/s)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x % 360:.0f}'))
    ax.set_title(f'HI4PI l-v diagram  {tag}  b = {glat_actual:.2f} deg')
    fig.colorbar(im, ax=ax, label='$T_B$ (K)')
    fig.tight_layout()
    out = safe_path(PLOT_DIR, f'HI4PI_lv_{tag}_b{glat_actual:+06.2f}_v{v_lo:+.0f}_{v_hi:+.0f}.png')
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return [out]


def _sky_map(values, header, wcs2d, title, cbar_label, fname, cmap, contours=False):
    glon_ax, glat_ax = world_axes(header, wcs2d)
    fig, ax = plt.subplots(figsize=(8, 7))
    im = ax.imshow(
        values, origin='lower', aspect='auto', cmap=cmap,
        extent=[glon_ax[0], glon_ax[-1], glat_ax[0], glat_ax[-1]])
    ax.set_xlabel('Galactic longitude $l$ (deg)')
    ax.set_ylabel('Galactic latitude $b$ (deg)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x % 360:.0f}'))
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    if contours:
        n_levels = 8
        vmin, vmax = np.nanmin(values), np.nanmax(values)
        levels = np.linspace(vmin, vmax, n_levels + 2)[1:-1]
        ax.contour(
            values, levels=levels, origin='lower', colors='white',
            linewidths=0.6, alpha=0.7,
            extent=[glon_ax[0], glon_ax[-1], glat_ax[0], glat_ax[-1]])
    fig.tight_layout()
    out = safe_path(PLOT_DIR, fname)
    fig.savefig(out, dpi=130)
    plt.close(fig)
    return out


def plot_moment0(data, header, wcs2d, v_kms, tag, vrange=None, contours=False):
    if vrange is None:
        vrange = ask_velocity_range(v_kms)
    v_lo, v_hi = vrange
    sel = (v_kms >= v_lo) & (v_kms <= v_hi)
    dv = abs(float(np.median(np.diff(v_kms))))
    mom0 = np.nansum(data[sel, :, :], axis=0) * dv   # K km/s
    nhi = mom0 * NHI_FACTOR
    out = _sky_map(
        nhi, header, wcs2d,
        f'HI4PI column density  {tag}  ({v_lo:.0f} to {v_hi:.0f} km/s)',
        r'$N_{HI}$ (cm$^{-2}$)',
        f'HI4PI_mom0_{tag}_v{v_lo:+.0f}_{v_hi:+.0f}.png', 'viridis', contours=contours)
    return [out]


def plot_peak_temp(data, header, wcs2d, v_kms, tag, contours=False):
    peak = np.nanmax(data, axis=0)
    out = _sky_map(
        peak, header, wcs2d, f'HI4PI peak brightness temperature  {tag}',
        r'$T_B^{peak}$ (K)', f'HI4PI_peakT_{tag}.png', 'inferno', contours=contours)
    return [out]


def plot_peak_velocity(data, header, wcs2d, v_kms, tag, contours=False):
    idx = np.nanargmax(data, axis=0)
    vpeak = v_kms[idx]
    out = _sky_map(
        vpeak, header, wcs2d, f'HI4PI velocity of peak emission  {tag}',
        'LSRK velocity of peak (km/s)', f'HI4PI_peakV_{tag}.png', 'coolwarm', contours=contours)
    return [out]


# ----------------------------------------------------------------- main ----

def main():
    print('=' * 62)
    print('  HI4PI visualisation suite - HALO project')
    print('=' * 62)

    fits_path = get_fits_path()
    if not os.path.isfile(fits_path):
        print(f'File not found: {fits_path}')
        sys.exit(1)

    tag = os.path.splitext(os.path.basename(fits_path))[0]

    global FRO_DIR, PLOT_DIR
    base_dir = os.path.join(SURVEYS_ROOT, f'HI4PI_{tag}')
    FRO_DIR = os.path.join(base_dir, 'FRO')
    PLOT_DIR = os.path.join(base_dir, 'Plots')
    os.makedirs(PLOT_DIR, exist_ok=True)

    print(f'\nLoading cube: {fits_path}')
    data, header = load_cube(fits_path)
    wcs2d = WCS(header).sub([1, 2])
    v_kms = velocity_axis_kms(header)
    print(f'Cube shape (vel, glat, glon): {data.shape}')
    print(f'Velocity range: {v_kms.min():.1f} .. {v_kms.max():.1f} km/s')
    print(f'Plots directory: {PLOT_DIR}')

    menu = """
  1) Spectra of the extracted .fro points
  2) l-v diagram at a chosen GLat
  3) Column density map (moment 0)
  4) Peak brightness temperature map
  5) Velocity-of-peak map
  6) All of the above
  0) Exit
"""
    while True:
        print(menu)
        choice = input('Choice: ').strip()
        contours = False
        if choice in ('3', '4', '5', '6'):
            contours = input('  Add contour lines? [y/N]: ').strip().lower() == 'y'
        written = []

        if choice == '0':
            break
        elif choice == '1':
            written = plot_extracted_spectra(tag)
        elif choice == '2':
            written = plot_lv(data, header, wcs2d, v_kms, tag)
        elif choice == '3':
            written = plot_moment0(data, header, wcs2d, v_kms, tag, contours=contours)
        elif choice == '4':
            written = plot_peak_temp(data, header, wcs2d, v_kms, tag, contours=contours)
        elif choice == '5':
            written = plot_peak_velocity(data, header, wcs2d, v_kms, tag, contours=contours)
        elif choice == '6':
            written += plot_extracted_spectra(tag)
            written += plot_lv(data, header, wcs2d, v_kms, tag)
            written += plot_moment0(data, header, wcs2d, v_kms, tag, contours=contours)
            written += plot_peak_temp(data, header, wcs2d, v_kms, tag, contours=contours)
            written += plot_peak_velocity(data, header, wcs2d, v_kms, tag, contours=contours)
        else:
            print('  Unknown choice.')
            continue

        for path in written:
            print(f'  written: {os.path.basename(path)}')
            try:
                subprocess.Popen(
                    ['xdg-open', path],
                    stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            except Exception as e:
                print(f'    (could not open viewer: {e})')

    print('\nDone.')


if __name__ == '__main__':
    main()