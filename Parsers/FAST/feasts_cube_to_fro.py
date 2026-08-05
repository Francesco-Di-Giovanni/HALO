#!/usr/bin/env python3
"""
feasts_cube_to_fro.py — FEASTS FAST HI cube to FRO format v17 parser.

Extracts spectra from FEASTS reduced HI data cubes (RA x Dec x VRAD)
at interactively chosen positions, writing one .fro file per session.

Input:  FEASTS FITS cube (e.g. NGC628-cube_cut.fits)
        Axes: NAXIS1=RA, NAXIS2=Dec, NAXIS3=VRAD [m/s]
        BUNIT: Jy/beam
        RESTFRQ: 1420405752.0 Hz

Output: ~/FRO/Parsers/FAST/Surveys/<galaxy>/<galaxy>_feasts.fro
        Plots in ~/FRO/Parsers/FAST/Surveys/<galaxy>/Plots/

Usage:  python3 feasts_cube_to_fro.py <path_to_cube.fits>
        or:     python3 feasts_cube_to_fro.py   (prompts for path)

Author: Francesco Di Giovanni (Claude AI assisted) - FRO/HALO Project
"""

import sys
import os
import subprocess
import numpy as np
from datetime import datetime, timezone
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Galactic
import astropy.units as u

# FRO format
sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
from fro_format_v17 import create_fro_v17_structure, append_raw_spectrum, UNKNOWN_EPOCH

# Physical constants
C_MS = 299792458.0   # speed of light [m/s]


def vrad_to_freq(v_ms, rest_freq_hz):
    """Convert radio velocity [m/s] to sky frequency [Hz]."""
    return rest_freq_hz * (1.0 - v_ms / C_MS)


def build_freq_axis(header):
    """Build full frequency axis [Hz] from VRAD WCS keywords."""
    naxis3 = header['NAXIS3']
    crpix3 = header['CRPIX3']
    crval3 = header['CRVAL3']   # m/s
    cdelt3 = header['CDELT3']   # m/s
    rest_freq = header['RESTFRQ']  # Hz
    v_array = crval3 + (np.arange(1, naxis3 + 1) - crpix3) * cdelt3
    return vrad_to_freq(v_array, rest_freq)


def radec_to_gal(ra_deg, dec_deg):
    """Convert RA/Dec [deg] to Galactic lon/lat [deg]."""
    c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
    g = c.galactic
    return g.l.deg, g.b.deg


def pixel_to_radec(wcs2d, ix, iy):
    """Convert 0-based pixel (ix, iy) to RA, Dec [deg]."""
    ra, dec = wcs2d.all_pix2world([[ix, iy]], 0)[0]
    return float(ra), float(dec)


def radec_to_pixel(wcs2d, ra_deg, dec_deg):
    """Convert RA, Dec [deg] to nearest 0-based pixel (ix, iy)."""
    px, py = wcs2d.all_world2pix([[ra_deg, dec_deg]], 0)[0]
    return int(round(px)), int(round(py))


def extract_spectrum(data, ix, iy):
    """Extract spectrum along VRAD axis at pixel (ix, iy). Returns 1D array."""
    return data[:, iy, ix].astype(np.float64)


def main():
    # --- Input file ---
    if len(sys.argv) > 1:
        cube_path = sys.argv[1]
    else:
        cube_path = input("Path to FEASTS cube FITS file: ").strip()
    cube_path = os.path.expanduser(cube_path)
    if not os.path.isfile(cube_path):
        print(f"ERROR: file not found: {cube_path}")
        sys.exit(1)

    # --- Load cube ---
    print(f"Loading {cube_path} ...")
    with fits.open(cube_path, memmap=True) as hdul:
        header = hdul[0].header
        data = hdul[0].data  # shape: (NAXIS3, NAXIS2, NAXIS1)

    galaxy = header.get('OBJECT', os.path.basename(cube_path).split('-')[0])
    telescope = header.get('TELESCOP', 'FAST')
    date_obs = header.get('DATE-OBS', 'unknown')
    bunit = header.get('BUNIT', 'Jy/beam')
    bmaj_deg = header.get('BMAJ', np.nan)
    rest_freq_hz = header.get('RESTFRQ', 1420405752.0)
    specsys = header.get('SPECSYS', 'BARYCENT')

    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']

    print(f"Galaxy : {galaxy}")
    print(f"Telescope: {telescope}")
    print(f"Date-Obs : {date_obs}")
    print(f"Cube shape: {data.shape}  (VRAD x Dec x RA)")
    print(f"BUNIT : {bunit}")
    print(f"Beam FWHM: {bmaj_deg * 60:.2f} arcmin")
    print(f"SPECSYS  : {specsys}")

    # --- Build WCS (2D spatial only) ---
    wcs_full = WCS(header)
    wcs2d = wcs_full.celestial

    # --- Build frequency axis ---
    freq_axis_hz = build_freq_axis(header)
    n_channels = len(freq_axis_hz)

    # --- Output directories ---
    surveys_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Surveys')
    galaxy_dir = os.path.join(surveys_dir, galaxy)
    fro_dir = os.path.join(galaxy_dir, 'FRO')
    plots_dir = os.path.join(galaxy_dir, 'Plots')
    for d in [fro_dir, plots_dir]:
        os.makedirs(d, exist_ok=True)

    fro_path = os.path.join(fro_dir, f'{galaxy}_feasts.fro')

    # --- Create FRO file ---
    obs_time = datetime.strptime(date_obs, '%Y-%m-%d').replace(tzinfo=timezone.utc)
    create_fro_v17_structure(
        filename=fro_path,
        frequency_axis_hz=freq_axis_hz,
        facility=telescope,
        instrument='19-beam L-band receiver',
        observation_id=f'FEASTS_{galaxy}_{date_obs}',
        survey='FEASTS',
        dataset=galaxy,
        provider='Wang et al. 2025 / kavli.pku.edu.cn/~jwang/FEASTS_data.html',
        observatory_name='FAST',
        latitude_deg=25.6530,
        longitude_deg=106.8566,
        elevation_m=1110.0,
        center_frequency_hz=rest_freq_hz,
        notes=(
            f'FEASTS reduced cube. Galaxy: {galaxy}. '
            f'BUNIT={bunit}. Beam={bmaj_deg * 60:.2f} arcmin. '
            f'SPECSYS={specsys}. RESTFRQ={rest_freq_hz:.0f} Hz. '
            f'Reduced by pipeline v1.2 (2025-07-31).'
        ),
    )
    print(f"\nFRO file created: {fro_path}")

    # --- Interactive extraction loop ---
    spectra_count = 0
    print("\nEnter RA/Dec [deg J2000] to extract spectra. Leave blank to finish.")
    while True:
        raw = input(f"\nRA [deg] (or blank to finish): ").strip()
        if not raw:
            break
        try:
            ra_req = float(raw)
            dec_req = float(input("Dec [deg]: ").strip())
        except ValueError:
            print("Invalid input, try again.")
            continue

        ix, iy = radec_to_pixel(wcs2d, ra_req, dec_req)

        # Boundary check
        if not (0 <= ix < naxis1 and 0 <= iy < naxis2):
            print(f"WARNING: pixel ({ix},{iy}) outside cube bounds ({naxis1}x{naxis2}). Skipping.")
            continue

        # Actual RA/Dec of extracted pixel
        ra_act, dec_act = pixel_to_radec(wcs2d, ix, iy)
        glon, glat = radec_to_gal(ra_act, dec_act)

        spectrum = extract_spectrum(data, ix, iy)

        nan_frac = np.sum(np.isnan(spectrum)) / n_channels
        if nan_frac > 0.5:
            print(f"WARNING: {nan_frac*100:.0f}% NaN channels at this position.")

        append_raw_spectrum(
            filename=fro_path,
            raw_spectrum=spectrum,
            utc_time=UNKNOWN_EPOCH,
            az_deg=None,
            elev_deg=None,
            glon_deg=glon,
            glat_deg=glat,
            ra_deg=ra_act,
            dec_deg=dec_act,
            integration_time_s=None,
            flag_track=1,
            flag_cal=0,
        )
        spectra_count += 1
        print(
            f"  Spectrum {spectra_count} extracted: "
            f"RA={ra_act:.4f} Dec={dec_act:.4f} "
            f"GLon={glon:.2f} GLat={glat:.2f}"
        )

    if spectra_count == 0:
        print("No spectra extracted. FRO file is empty.")
        sys.exit(0)

    print(f"\nTotal spectra written: {spectra_count}")
    print(f"FRO file: {fro_path}")

    # --- Offer to plot ---
    plot_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'plot_feasts.py')
    if os.path.isfile(plot_script):
        ans = input("\nGenerate plots now? [y/N]: ").strip().lower()
        if ans == 'y':
            subprocess.run([sys.executable, plot_script, cube_path], check=False)


if __name__ == '__main__':
    main()