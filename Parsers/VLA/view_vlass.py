# view_vlass.py
# Viewer for VLASS continuum FITS cutouts
# Produces jet colormap image with contours, WCS axes in sexagesimal RA/Dec
# Nothing hardcoded: file path from sys.argv or input()
#
# Author: Francesco Di Giovanni / HALO project
# Date: 2026-07-28

import os
import sys
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import AsinhStretch, ImageNormalize
from astropy.coordinates import SkyCoord
import astropy.units as u



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
SURVEYS_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Surveys')


def get_fits_path():
    if len(sys.argv) > 1:
        return os.path.expanduser(sys.argv[1])
    return os.path.expanduser(input('Path to VLASS FITS cutout: ').strip())


def load_image(filepath):
    with fits.open(filepath) as hdul:
        data = hdul[0].data.copy()
        header = hdul[0].header.copy()
    # Squeeze degenerate axes (Stokes, Freq) to get 2D
    data = np.squeeze(data)
    # Build 2D WCS from the celestial axes only
    wcs_full = WCS(header)
    wcs2d = wcs_full.celestial
    return data, header, wcs2d


def plot_continuum(data, header, wcs2d, fits_path, contours=True):
    tag = os.path.splitext(os.path.basename(fits_path))[0]
    obj = header.get('OBJECT', tag)
    freq_hz = header.get('CRVAL3', header.get('RESTFRQ', np.nan))
    freq_ghz = freq_hz / 1e9 if np.isfinite(freq_hz) else 0
    bunit = header.get('BUNIT', 'Jy/beam')
    bmaj_arcsec = header.get('BMAJ', 0) * 3600
    bmin_arcsec = header.get('BMIN', 0) * 3600
    bpa = header.get('BPA', 0)

    # Mask NaN
    data = np.where(np.isfinite(data), data, np.nan)

    # Normalization: arcsinh stretch like the CIRADA service
    vmin = np.nanpercentile(data, 0.5)
    vmax = np.nanpercentile(data, 99.9)
    norm = ImageNormalize(vmin=vmin, vmax=vmax, stretch=AsinhStretch(0.1))

    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection=wcs2d)

    im = ax.imshow(data, origin='lower', cmap='jet', norm=norm)

    # Contours
    if contours:
        peak = np.nanmax(data)
        rms = np.nanstd(data[data < np.nanpercentile(data, 50)])
        if rms > 0 and peak > 5 * rms:
            levels = rms * np.array([3, 5, 10, 20, 50, 100, 200, 500])
            levels = levels[levels < peak]
            if len(levels) > 0:
                ax.contour(data, levels=levels, colors='white',
                           linewidths=0.5, alpha=0.8, transform=ax.get_transform(wcs2d))

    # Colorbar
    cbar = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.02)
    cbar.set_label(bunit, fontsize=11)

    # Title with survey info
    title = f'VLASS  {obj}'
    if freq_ghz > 0:
        title += f'    {freq_ghz:.1f} GHz'
    ax.set_title(title, fontsize=13, pad=12)

    # Beam ellipse in lower-left corner
    if bmaj_arcsec > 0 and bmin_arcsec > 0:
        from matplotlib.patches import Ellipse
        pixel_scale = abs(header.get('CDELT1', 1)) * 3600  # arcsec/pixel
        beam_x = bmaj_arcsec / pixel_scale
        beam_y = bmin_arcsec / pixel_scale
        beam_patch = Ellipse(
            (15, 15), beam_x, beam_y, angle=bpa,
            edgecolor='white', facecolor='gray', alpha=0.7, lw=1)
        ax.add_patch(beam_patch)

    # White text/ticks on dark background
    ax.tick_params(colors='white', which='both')
    ax.coords[0].set_ticks(spacing=30 * u.arcsec)
    ax.coords[1].set_ticks(spacing=30 * u.arcsec)
    ax.coords[0].set_major_formatter("hh:mm:ss")
    ax.coords[1].set_major_formatter("dd:mm:ss")
    ax.coords[0].set_ticklabel(color='white')
    ax.coords[1].set_ticklabel(color='white')
    ax.coords[0].set_axislabel('RA (J2000)', color='white', fontsize=11)
    ax.coords[1].set_axislabel('Dec (J2000)', color='white', fontsize=11)
    ax.set_title(title, fontsize=13, pad=12, color='white')
    cbar.ax.yaxis.set_tick_params(color='white')
    cbar.ax.yaxis.label.set_color('white')
    plt.setp(cbar.ax.yaxis.get_ticklabels(), color='white')
    for spine in ax.spines.values():
        spine.set_edgecolor('white')

    fig.tight_layout()

    # Save
    base_dir = os.path.join(SURVEYS_ROOT, f'VLASS_{tag}')
    plot_dir = os.path.join(base_dir, 'Plots')
    os.makedirs(plot_dir, exist_ok=True)
    out_path = safe_path(plot_dir, f'VLASS_{tag}_continuum.png')
    fig.savefig(out_path, dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    return out_path


def main():
    print('=' * 60)
    print('  VLASS Continuum Viewer - HALO project')
    print('=' * 60)

    fits_path = get_fits_path()
    if not os.path.isfile(fits_path):
        print(f'File not found: {fits_path}')
        sys.exit(1)

    print(f'\nLoading: {fits_path}')
    data, header, wcs2d = load_image(fits_path)
    print(f'Image shape: {data.shape}')
    print(f'Object: {header.get("OBJECT", "unknown")}')
    print(f'Frequency: {header.get("CRVAL3", 0)/1e9:.2f} GHz')
    print(f'BUNIT: {header.get("BUNIT", "unknown")}')
    print(f'Beam: {header.get("BMAJ",0)*3600:.2f}" x {header.get("BMIN",0)*3600:.2f}"')

    contours = input('\nAdd contour lines? [Y/n]: ').strip().lower() != 'n'

    out_path = plot_continuum(data, header, wcs2d, fits_path, contours=contours)
    print(f'\nWritten: {out_path}')

    try:
        subprocess.Popen(['xdg-open', out_path],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    except Exception:
        pass

    print('Done.')


if __name__ == '__main__':
    main()