# hi4pi_to_fro.py
# Converts HI4PI survey FITS cubes (position-position-velocity) to FRO HDF5 format (v17)
# HI4PI = combination of EBHIS (Effelsberg 100m) and GASS (Parkes 64m) all-sky HI surveys
# Reference: HI4PI Collaboration 2016, A&A 594, A116
#
# Unlike single-dish parsers, this reads one spectrum at a time from a 3D cube,
# at galactic coordinates chosen interactively by the user (no fixed grid,
# no assumption about which points will ever be extracted).
#
# Facility/observatory assigned dynamically per point, based on declination
# (Dec > -5 deg -> Effelsberg/EBHIS component, Dec <= -5 deg -> Parkes/GASS component),
# since HI4PI is a combined all-sky product with no single physical site.
#
# Author: Francesco Di Giovanni / HALO project
# Date: 2026-07-27

import os
import sys
import subprocess
import numpy as np
from datetime import datetime, timezone
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord, Galactic, FK5
import astropy.units as u

sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
import fro_format_v17 as fro

C_MPS = 299792458.0  # m/s

# Effelsberg 100m coordinates
EFFELSBERG_LAT = 50.5248
EFFELSBERG_LON = 6.8836
EFFELSBERG_ALT = 369.0

# Parkes 64m coordinates
PARKES_LAT = -32.9984
PARKES_LON = 148.2621
PARKES_ALT = 392.0

DEC_SPLIT_DEG = -5.0  # EBHIS/GASS survey boundary in declination


def get_fits_path():
    """Return the path to the HI4PI FITS cube: from argv[1] if given,
    otherwise ask interactively. Never hardcoded to a specific tile."""
    if len(sys.argv) > 1:
        return sys.argv[1]
    return input('Path to HI4PI FITS cube: ').strip()


def load_cube(filepath):
    with fits.open(filepath) as hdul:
        data = hdul[0].data.astype(np.float32)
        header = hdul[0].header.copy()
    return data, header


def build_velocity_axis(header):
    """Reconstruct the LSRK radio velocity axis (m/s) from the FITS header.
    Nothing hardcoded: CRPIX3/CRVAL3/CDELT3 read dynamically."""
    n = header['NAXIS3']
    crpix3 = header['CRPIX3']
    crval3 = header['CRVAL3']
    cdelt3 = header['CDELT3']
    pix = np.arange(1, n + 1)
    v_mps = crval3 + (pix - crpix3) * cdelt3  # CUNIT3 = m/s
    return v_mps


def velocity_to_freq(v_mps, restfreq_hz):
    return restfreq_hz * (1.0 - v_mps / C_MPS)


def extract_spectrum(data, header, wcs2d, glon_req, glat_req):
    """Extract the spectrum at the pixel nearest to (glon_req, glat_req).
    Returns the spectrum plus the ACTUAL pixel-center coordinates,
    for honest metadata (requested vs actual)."""
    ix, iy = wcs2d.wcs_world2pix(glon_req, glat_req, 0)
    ix = int(round(float(ix)))
    iy = int(round(float(iy)))

    nx = header['NAXIS1']
    ny = header['NAXIS2']
    if not (0 <= ix < nx and 0 <= iy < ny):
        raise ValueError(
            f'Requested (GLon={glon_req}, GLat={glat_req}) falls outside '
            f'the cube coverage (pixel {ix},{iy} out of {nx}x{ny}).'
        )

    glon_actual, glat_actual = wcs2d.wcs_pix2world(ix, iy, 0)
    glon_actual = float(glon_actual)
    glat_actual = float(glat_actual)

    spectrum = data[:, iy, ix].astype(np.float32)
    return spectrum, ix, iy, glon_actual, glat_actual


def convert_point(data, header, wcs2d, glon_req, glat_req, output_dir, tile_name):
    restfreq_hz = header['RESTFRQ']

    spectrum, ix, iy, glon_actual, glat_actual = extract_spectrum(
        data, header, wcs2d, glon_req, glat_req)

    v_mps = build_velocity_axis(header)
    freq_hz = velocity_to_freq(v_mps, restfreq_hz)

    sort_idx = np.argsort(freq_hz)
    freq_hz_sorted = freq_hz[sort_idx]
    spectrum_sorted = spectrum[sort_idx]

    gc = SkyCoord(l=glon_actual * u.deg, b=glat_actual * u.deg, frame=Galactic())
    eq = gc.transform_to(FK5(equinox='J2000'))
    ra_deg = eq.ra.deg
    dec_deg = eq.dec.deg

    if dec_deg > DEC_SPLIT_DEG:
        facility = 'Effelsberg 100m'
        instrument = 'EBHIS component of HI4PI'
        obs_lat, obs_lon, obs_alt = EFFELSBERG_LAT, EFFELSBERG_LON, EFFELSBERG_ALT
    else:
        facility = 'Parkes 64m'
        instrument = 'GASS component of HI4PI'
        obs_lat, obs_lon, obs_alt = PARKES_LAT, PARKES_LON, PARKES_ALT

    out_name = f'HI4PI_G{glon_actual:06.2f}{glat_actual:+06.2f}_FRO.fro'
    out_path = os.path.join(output_dir, out_name)

    notes = (
        f'HI4PI all-sky HI survey (HI4PI Collaboration 2016, A&A 594, A116). '
        f'Source tile: {tile_name}. '
        f'Requested coordinates: GLon={glon_req:.4f}, GLat={glat_req:.4f}. '
        f'Actual extracted pixel-center coordinates: GLon={glon_actual:.4f}, '
        f'GLat={glat_actual:.4f} (pixel {ix},{iy}; grid step '
        f'{abs(header["CDELT1"]):.5f} deg). '
        f'Facility assigned dynamically by declination (Dec={dec_deg:.3f} deg, '
        f'split at Dec={DEC_SPLIT_DEG} deg): {facility}. '
        f'Frequency axis reconstructed from CRVAL3/CRPIX3/CDELT3 header keywords '
        f'(LSRK radio velocity, m/s) using f = f0*(1-v/c), '
        f'f0 = RESTFRQ = {restfreq_hz:.2f} Hz (from FITS header). '
        f'Observation epoch: not defined for this product. Each cube pixel '
        f'combines multiple original observations from the two component '
        f'surveys (EBHIS, Effelsberg, observed ~2008-2013; GASS, Parkes, '
        f'~2005-2006), resampled onto a common grid; Time/utc is therefore '
        f'set to the explicit UNKNOWN_EPOCH sentinel rather than to any date.'
    )

    fro.create_fro_v17_structure(
        filename=out_path,
        frequency_axis_hz=freq_hz_sorted,
        facility=facility,
        instrument=instrument,
        observation_id=f'{tile_name}_px{ix}_{iy}',
        survey='HI4PI',
        dataset=f'HI4PI GLon{glon_actual:.2f} GLat{glat_actual:.2f}',
        provider='HI4PI Collaboration (AIfA/MPIfR Bonn; ATNF Sydney)',
        import_date=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        observatory_name=facility,
        latitude_deg=obs_lat,
        longitude_deg=obs_lon,
        elevation_m=obs_alt,
        receiver='L-band',
        center_frequency_hz=restfreq_hz,
        notes=notes,
        pointing_source='computed',
        temperature_calibration_note=(
            'Brightness temperature (K) on the HI4PI calibrated scale '
            '(BUNIT=K in source FITS; Effelsberg/Parkes absolute calibration). '
            'This is a calibrated product, not an assumed-Tsys nominal scale.'
        ),
    )

    fro.append_raw_spectrum(
        filename=out_path,
        raw_spectrum=spectrum_sorted,
        utc_time=fro.UNKNOWN_EPOCH,
        az_deg=float('nan'),
        elev_deg=float('nan'),
        ra_deg=ra_deg,
        dec_deg=dec_deg,
        glon_deg=glon_actual,
        glat_deg=glat_actual,
        integration_time_s=float('nan'),
        polarization='I',
        flag_track=-1,
        flag_cal=0,
    )

    return out_path, dec_deg, facility


def main():
    print('=' * 60)
    print('  HI4PI to FRO converter - HALO project')
    print('  Interactive point extraction from position-position-velocity cube')
    print('=' * 60)

    fits_path = os.path.expanduser(get_fits_path())
    if not os.path.isfile(fits_path):
        print(f'File not found: {fits_path}')
        sys.exit(1)

    tile_name = os.path.splitext(os.path.basename(fits_path))[0]

    print(f'\nLoading cube: {fits_path}')
    data, header = load_cube(fits_path)
    wcs2d = WCS(header).sub([1, 2])

    print(f'Cube shape (vel, glat, glon): {data.shape}')
    print(f'RESTFRQ: {header["RESTFRQ"]:.2f} Hz')
    print(f'Grid step: {abs(header["CDELT1"]):.5f} deg')

    # Working directory derived from the cube file itself: nothing fixed.
    output_dir = os.path.expanduser(
        f'~/FRO/Parsers/HI4PI/Surveys/HI4PI_{tile_name}/FRO')
    os.makedirs(output_dir, exist_ok=True)
    print(f'Output directory: {output_dir}\n')

    while True:
        try:
            glon_req = float(input('GLon (deg): '))
            glat_req = float(input('GLat (deg): '))
        except ValueError:
            print('Invalid input, try again.\n')
            continue

        try:
            out_path, dec_deg, facility = convert_point(
                data, header, wcs2d, glon_req, glat_req, output_dir, tile_name)
            print(f'  OK -> {os.path.basename(out_path)}  '
                  f'(Dec={dec_deg:.2f}, {facility})')
        except Exception as e:
            print(f'  ERROR: {e}')

        again = input('\nExtract another point? [y/N]: ').strip().lower()
        if again != 'y':
            break

    print('\nDone.')

    # Hand over to the plotting suite, on the same cube, if wanted.
    plotter = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'plot_hi4pi.py')
    if os.path.isfile(plotter):
        if input('\nGenerate plots now? [y/N]: ').strip().lower() == 'y':
            subprocess.call([sys.executable, plotter, fits_path])


if __name__ == '__main__':
    main()
