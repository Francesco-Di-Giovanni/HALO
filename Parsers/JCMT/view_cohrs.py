# view_cohrs.py
# Visualisation suite for COHRS 12CO(3-2) position-position-velocity cubes
# (JCMT / HARP, CO High-Resolution Survey of the Galactic Plane) and for the
# .fro spectra extracted from them by cohrs_to_fro.py.
#
# Nothing is hardcoded: cube path, latitude of the l-v cut and velocity
# integration range are all supplied at runtime. The working directory is
# derived automatically from the cube file name.
#
# All figures are written to Surveys/COHRS_<tag>/Plots/
#
# Author: Francesco Di Giovanni / HALO project
# Date: 2026-07-31

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
import warnings
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning

C_MPS = 299792458.0

SURVEYS_ROOT = os.path.expanduser('~/FRO/Parsers/JCMT/Surveys')

# Set at runtime from the chosen cube file: no fixed working directory.
FRO_DIR = ''
PLOT_DIR = ''
PLOT_DIR_CONTOUR = ''
PLOT_TIMESTAMP = ''

def get_plot_dir(contours=False):
    """Return the correct output directory based on contour flag."""
    return PLOT_DIR_CONTOUR if contours else PLOT_DIR

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
    surveys_dir = os.path.expanduser('~/FRO/Parsers/JCMT/Surveys/')
    result = subprocess.run(
        ['zenity', '--file-selection',
         '--title=Select COHRS FITS cube',
         '--file-filter=FITS files (*.fit *.fits)|*.fit *.fits',
         f'--filename={surveys_dir}'],
        capture_output=True, text=True)
    fits_path = result.stdout.strip()
    if not fits_path:
        print('No file selected, exiting.')
        sys.exit(1)
    return os.path.expanduser(fits_path)


def load_cube(filepath):
    """Load the PRIMARY data cube and header. Byte order is normalised to
    native so NumPy 2.0 does not choke on the big-endian FITS array."""
    with fits.open(filepath) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
        header = hdul[0].header.copy()
    return data, header


def velocity_axis_kms(header):
    """LSRK radio velocity axis in km/s, taken from the header only.
    COHRS stores CTYPE3=VRAD with CUNIT3=km/s, so the axis is already a
    velocity: no frequency conversion is needed. If a future cube uses m/s
    the CUNIT3 keyword is honoured."""
    n = header['NAXIS3']
    pix = np.arange(1, n + 1)
    v = header['CRVAL3'] + (pix - header['CRPIX3']) * header['CDELT3']
    cunit = str(header.get('CUNIT3', 'km/s')).strip().lower()
    if cunit in ('m/s', 'ms-1', 'm s-1'):
        v = v / 1000.0
    return v


def rest_frequency_hz(header):
    for k in ('RESTFRQ', 'RESTFREQ'):
        if k in header:
            return float(header[k])
    return np.nan


def world_axes(header, wcs2d):
    """GLon and GLat values along the two spatial axes, read from the WCS.
    Handles the 360/0 degree wrap for tiles crossing l=0 (not the case for
    the G30 tile, but kept for generality)."""
    nx = header['NAXIS1']
    ny = header['NAXIS2']
    ix = np.arange(nx)
    iy = np.arange(ny)
    glon, _ = wcs2d.wcs_pix2world(ix, np.zeros(nx), 0)
    _, glat = wcs2d.wcs_pix2world(np.zeros(ny), iy, 0)
    glon = np.asarray(glon, dtype=float)
    glat = np.asarray(glat, dtype=float)
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


def ask_colormap(default='inferno'):
    txt = input(f'  Colormap [inferno/viridis/jet/coolwarm, blank={default}]: ').strip()
    return txt if txt else default


def clip_limits(values, lo_pct=1.0, hi_pct=99.5):
    """Robust colour limits from percentiles of the finite values, so a few
    very bright pixels (e.g. W43) do not wash out the whole map."""
    fin = values[np.isfinite(values)]
    if fin.size == 0:
        return None, None
    return float(np.percentile(fin, lo_pct)), float(np.percentile(fin, hi_pct))


# ---------------------------------------------------------------- plots ----

def plot_extracted_spectra(tag, contours=False):
    """One figure per .fro previously produced by cohrs_to_fro.py."""
    files = sorted(glob.glob(os.path.join(FRO_DIR, '*.fro')))
    if not files:
        print(f'  No .fro files found in {FRO_DIR}')
        print('  (run cohrs_to_fro.py first to extract spectra)')
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
            restfreq = 345795989900.0
        v_kms = C_MPS * (1.0 - freq / restfreq) / 1000.0
        order = np.argsort(v_kms)

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.plot(v_kms[order], spec[order], lw=1.0, color='#8a1f9c')
        ax.set_xlabel('LSRK velocity (km/s)')
        restfreq_ghz = restfreq / 1e9
        secax = ax.secondary_xaxis('top', functions=(lambda v: restfreq_ghz * (1 - v / 299792.458), lambda f: (1 - f / restfreq_ghz) * 299792.458))
        secax.set_xlabel('Frequency (GHz)')
        ax.set_ylabel(r'Antenna temperature $T_A^*$ (K)')
        ax.set_title(f'COHRS 12CO(3-2) spectrum  GLon={glon:.3f}  GLat={glat:.3f}')
        ax.grid(alpha=0.3)
        from scipy.signal import find_peaks
        rms = np.nanstd(spec[np.abs(v_kms) > 150])
        all_peaks, _ = find_peaks(spec, prominence=2*rms, distance=3)
        all_peaks = all_peaks[spec[all_peaks] > 4.0]
        top_n = 3
        if len(all_peaks) > top_n:
            all_peaks = all_peaks[np.argsort(spec[all_peaks])[::-1][:top_n]]
        peaks = all_peaks[np.argsort(v_kms[all_peaks])]
        v_span = v_kms.max() - v_kms.min()
        y_max = np.nanmax(spec)
        ax.set_ylim(top=y_max * 1.45)
        for ii, pi in enumerate(peaks):
            vp = v_kms[pi]
            fp = restfreq_ghz * (1 - vp / 299792.458)
            label = f"v={vp:.1f} km/s" + chr(10) + f"{fp:.4f} GHz"
            side = -1 if ii == 0 else 1
            x_off = side * v_span * 0.10
            y_top = spec[pi] + (y_max - spec.min()) * 0.25
            ax.annotate(label, xy=(vp, spec[pi]), xytext=(vp + x_off, y_top), fontsize=7, color="red", ha="center", va="bottom", arrowprops=dict(arrowstyle="->", color="red", lw=0.8))
            ax.plot(vp, spec[pi], marker="v", color="red", markersize=6)
        fig.tight_layout()
        fname_sp = f'COHRS_spectrum_G{glon:07.3f}{glat:+07.3f}_{PLOT_TIMESTAMP}{"_contour" if contours else ""}.svg'
        out = safe_path(get_plot_dir(contours), fname_sp)
        fig.savefig(out, format="svg")
        plt.close(fig)
        written.append(out)
    return written


def plot_lv(data, header, wcs2d, v_kms, tag, contours=False, glat_txt=None, vrange=None):
    glon_ax, glat_ax = world_axes(header, wcs2d)
    if glat_txt is None:
        print(f'  GLat available: {glat_ax.min():.3f} .. {glat_ax.max():.3f}')
        txt = input('  GLat of the l-v cut (deg, blank = integrate over all b): ').strip()
    else:
        txt = glat_txt
    if vrange is None:
        v_lo, v_hi = ask_velocity_range(v_kms)
    else:
        v_lo, v_hi = vrange
    sel = (v_kms >= v_lo) & (v_kms <= v_hi)
    v_cut = v_kms[sel]

    if txt == '':
        # collapse (average) over latitude -> classic COHRS l-v map
        slab = np.nanmean(data[sel, :, :], axis=1)  # (velocity, glon)
        lat_label = 'b integrated'
        lat_tag = 'ball'
    else:
        try:
            glat_req = float(txt)
        except ValueError:
            print('  Invalid latitude, skipping.')
            return []
        iy = int(np.argmin(np.abs(glat_ax - glat_req)))
        glat_actual = glat_ax[iy]
        slab = data[sel, iy, :]  # (velocity, glon)
        lat_label = f'b = {glat_actual:.3f} deg'
        lat_tag = f'b{glat_actual:+06.2f}'

    vmin, vmax = clip_limits(slab)
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(
        slab, origin='lower', aspect='auto', cmap='inferno',
        vmin=vmin, vmax=vmax,
        extent=[glon_ax[0], glon_ax[-1], v_cut[0], v_cut[-1]])
    ax.set_xlabel('Galactic longitude $l$ (deg)')
    ax.set_ylabel('LSRK velocity (km/s)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x % 360:.1f}'))
    ax.set_title(f'COHRS 12CO(3-2) l-v diagram  {tag}  ({lat_label})')
    fig.colorbar(im, ax=ax, label=r'$T_A^*$ (K)')
    fig.tight_layout()
    if contours:
        levels = np.linspace(np.nanpercentile(slab, 10),
                             np.nanpercentile(slab, 99), 8)
        ax.contour(slab, origin='lower', levels=levels, colors='white',
                   linewidths=0.5, alpha=0.7,
                   extent=[glon_ax[0], glon_ax[-1], v_cut[0], v_cut[-1]])
    fig.tight_layout()
    fname_lv = f'lv_{tag}_{lat_tag}_v{v_lo:+.0f}_{v_hi:+.0f}_{PLOT_TIMESTAMP}{"_contour" if contours else ""}.png'
    out = safe_path(get_plot_dir(contours), fname_lv)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return [out]


def _sky_map(values, header, wcs2d, title, cbar_label, fname, cmap,
             contours=False, clip=True):
    glon_ax, glat_ax = world_axes(header, wcs2d)
    if clip:
        vmin, vmax = clip_limits(values)
    else:
        vmin = vmax = None
    fig, ax = plt.subplots(figsize=(11, 5))
    im = ax.imshow(
        values, origin='lower', aspect='auto', cmap=cmap,
        vmin=vmin, vmax=vmax, interpolation='bilinear',
        extent=[glon_ax[0], glon_ax[-1], glat_ax[0], glat_ax[-1]])
    ax.set_xlabel('Galactic longitude $l$ (deg)')
    ax.set_ylabel('Galactic latitude $b$ (deg)')
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x % 360:.1f}'))
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label=cbar_label)
    if contours:
        n_levels = 8
        fin = values[np.isfinite(values)]
        if fin.size:
            lo, hi = np.percentile(fin, [50, 99])
            levels = np.linspace(lo, hi, n_levels)
            ax.contour(
                values, levels=levels, origin='lower', colors='white',
                linewidths=0.5, alpha=0.7,
                extent=[glon_ax[0], glon_ax[-1], glat_ax[0], glat_ax[-1]])
    fig.tight_layout()
    out = safe_path(get_plot_dir(contours), fname)
    fig.savefig(out, dpi=300)
    plt.close(fig)
    return out


def plot_moment0(data, header, wcs2d, v_kms, tag, vrange=None,
                 contours=False, cmap='inferno'):
    """Velocity-integrated 12CO(3-2) intensity W(CO) in K km/s.
    We deliberately do NOT convert to a column density: that would require
    an assumed CO-to-H2 conversion factor (X_CO), which is a strong and
    uncertain assumption. W(CO) is the honest, directly measured quantity."""
    if vrange is None:
        vrange = ask_velocity_range(v_kms)
    v_lo, v_hi = vrange
    sel = (v_kms >= v_lo) & (v_kms <= v_hi)
    dv = abs(float(np.median(np.diff(v_kms))))
    wco = np.nansum(data[sel, :, :], axis=0) * dv   # K km/s
    out = _sky_map(
        wco, header, wcs2d,
        f'COHRS 12CO(3-2) integrated intensity  {tag}  '
        f'({v_lo:.0f} to {v_hi:.0f} km/s)',
        r'$W_{CO}$ (K km/s)',
        f'mom0_{tag}_v{v_lo:+.0f}_{v_hi:+.0f}_{PLOT_TIMESTAMP}{"_contour" if contours else ""}.png', 'inferno',
        contours=contours)
    return [out]


def plot_peak_temp(data, header, wcs2d, v_kms, tag, contours=False, cmap='inferno'):
    peak = np.nanmax(data, axis=0)
    out = _sky_map(
        peak, header, wcs2d,
        f'COHRS 12CO(3-2) peak antenna temperature  {tag}',
        r'$T_A^{*,peak}$ (K)', f'peakT_{tag}_{PLOT_TIMESTAMP}{"_contour" if contours else ""}.png', 'inferno',
        contours=contours)
    return [out]


def plot_peak_velocity(data, header, wcs2d, v_kms, tag, contours=False,
                       cmap='coolwarm'):
    """Velocity of the peak channel. Pixels whose peak is below a small
    threshold are masked to NaN, otherwise noise-dominated lines of sight
    paint random velocities across the map."""
    peak = np.nanmax(data, axis=0)
    idx = np.nanargmax(data, axis=0)
    vpeak = v_kms[idx].astype(np.float32)
    thr = float(np.nanpercentile(peak[np.isfinite(peak)], 30))
    vpeak[peak < thr] = np.nan
    out = _sky_map(
        vpeak, header, wcs2d,
        f'COHRS 12CO(3-2) velocity of peak emission  {tag}',
        'LSRK velocity of peak (km/s)', f'peakV_{tag}_{PLOT_TIMESTAMP}{"_contour" if contours else ""}.png', cmap,
        contours=contours, clip=False)
    return [out]


# ----------------------------------------------------------------- main ----

def main():
    print('=' * 62)
    print('  COHRS 12CO(3-2) visualisation suite - HALO project')
    print('=' * 62)

    fits_path = get_fits_path()
    if not os.path.isfile(fits_path):
        print(f'File not found: {fits_path}')
        sys.exit(1)

    tag = os.path.splitext(os.path.basename(fits_path))[0]

    global FRO_DIR, PLOT_DIR, PLOT_DIR_CONTOUR
    if len(sys.argv) > 2:
        epoch_dir = os.path.expanduser(sys.argv[2])
    else:
        base_dir = os.path.join(SURVEYS_ROOT, tag)
        result = subprocess.run(
            ['zenity', '--file-selection', '--directory',
             '--title=Select epoch directory',
             f'--filename={base_dir}/'],
            capture_output=True, text=True)
        epoch_dir = result.stdout.strip()
        if not epoch_dir:
            print('No epoch directory selected, exiting.')
            sys.exit(1)
    FRO_DIR = os.path.join(epoch_dir, 'FRO')
    PLOT_DIR = os.path.join(epoch_dir, 'Plots')
    PLOT_DIR_CONTOUR = os.path.join(epoch_dir, 'Plots', 'Plots_contour')
    os.makedirs(PLOT_DIR, exist_ok=True)
    os.makedirs(PLOT_DIR_CONTOUR, exist_ok=True)
    # Extract timestamp from epoch directory name (last 15 chars: YYYYMMDD_HHMMSS)
    global PLOT_TIMESTAMP
    PLOT_TIMESTAMP = os.path.basename(epoch_dir)[-15:]

    print(f'\nLoading cube: {fits_path}')
    data, header = load_cube(fits_path)
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always', FITSFixedWarning)
        wcs2d = WCS(header).sub([1, 2])
    if caught:
        fixes = ', '.join(sorted(set(
            w.message.args[0].split("'")[1] for w in caught
            if issubclass(w.category, FITSFixedWarning))))
        print(f'  Header corrections applied by astropy WCS ({fixes}) '
              f'because the FITS header uses non-standard conventions - OK')
    v_kms = velocity_axis_kms(header)
    restfreq = rest_frequency_hz(header)
    print(f'Cube shape (vel, glat, glon): {data.shape}')
    print(f'Velocity range: {v_kms.min():.1f} .. {v_kms.max():.1f} km/s '
          f'({header["CDELT3"]:.3f} km/s/channel)')
    if np.isfinite(restfreq):
        print(f'Rest frequency: {restfreq/1e9:.6f} GHz  '
              f'(line: {header.get("MOLECULE", "12CO(3-2)")})')
    print(f'Telescope/instrument: {header.get("TELESCOP","?")}/'
          f'{header.get("INSTRUME","?")}')
    print(f'Plots directory: {PLOT_DIR}')

    menu = """
  1) Spectra of the extracted .fro points
  2) l-v diagram (chosen GLat, or integrated over b)
  3) Integrated intensity map W(CO) (moment 0)
  4) Peak antenna temperature map
  5) Velocity-of-peak map
  6) All of the above
  0) Exit
"""
    while True:
        print(menu)
        choice = input('Choice: ').strip()
        cmap = None
        contour_choice = 'n'
        if choice in ('1', '2', '3', '4', '5', '6'):
            print('  Contour lines?')
            print('    1) No contours')
            print('    2) With contours')
            print('    3) Both')
            contour_choice = input('  Choice [1]: ').strip()
        contour_modes = []
        if contour_choice == '2':
            contour_modes = [True]
        elif contour_choice == '3':
            contour_modes = [False, True]
        else:
            contour_modes = [False]
        if choice in ('3', '4', '5'):
            cmap = ask_colormap({'3': 'viridis', '4': 'inferno',
                                 '5': 'coolwarm'}[choice])
        written = []

        if choice == '0':
            break
        elif choice not in ('1','2','3','4','5','6'):
            print('  Unknown choice.')
            continue
        if choice in ('2', '6'):
            glon_ax, glat_ax = world_axes(header, wcs2d)
            print(f'  GLat available: {glat_ax.min():.3f} .. {glat_ax.max():.3f}')
            glat_txt = input('  GLat of the l-v cut (deg, blank = integrate over all b): ').strip()
            lv_vrange = ask_velocity_range(v_kms)
        if choice == '6':
            vrange = lv_vrange
        for contours in contour_modes:
            if choice == '1':
                written += plot_extracted_spectra(tag, contours=contours)
            elif choice == '2':
                written += plot_lv(data, header, wcs2d, v_kms, tag, contours=contours,
                                   glat_txt=glat_txt, vrange=lv_vrange)
            elif choice == '3':
                written += plot_moment0(data, header, wcs2d, v_kms, tag,
                                       contours=contours, cmap=cmap)
            elif choice == '4':
                written += plot_peak_temp(data, header, wcs2d, v_kms, tag,
                                         contours=contours, cmap=cmap)
            elif choice == '5':
                written += plot_peak_velocity(data, header, wcs2d, v_kms, tag,
                                             contours=contours, cmap=cmap)
            elif choice == '6':
                written += plot_extracted_spectra(tag, contours=contours)
                written += plot_lv(data, header, wcs2d, v_kms, tag, contours=contours,
                                   glat_txt=glat_txt, vrange=lv_vrange)
                written += plot_moment0(data, header, wcs2d, v_kms, tag,
                                        vrange=vrange, contours=contours,
                                        cmap='inferno')
                written += plot_peak_temp(data, header, wcs2d, v_kms, tag,
                                          contours=contours, cmap='inferno')
                written += plot_peak_velocity(data, header, wcs2d, v_kms, tag,
                                              contours=contours, cmap='coolwarm')

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