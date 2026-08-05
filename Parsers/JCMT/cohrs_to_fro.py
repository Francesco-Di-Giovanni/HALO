# cohrs_to_fro.py
# Converts COHRS survey FITS cubes (position-position-velocity) to FRO HDF5 format (v17)
# COHRS = 12CO(3-2) High-Resolution Survey of the Galactic Plane (JCMT / HARP)
# Reference: Park et al. 2023, ApJS 264, 16 (Complete Data Release, DR2);
#            Dempsey, Thomas & Currie 2013, ApJS 209, 8 (DR1)
#
# Like hi4pi_to_fro.py, this reads one spectrum at a time from a 3D cube,
# at galactic coordinates chosen interactively by the user (no fixed grid,
# no assumption about which points will ever be extracted).
#
# Unlike HI4PI, the facility is fixed (JCMT / HARP): COHRS is a single-telescope
# survey. The spectral axis of the cube is already a radial velocity (CTYPE3=VRAD,
# CUNIT3=km/s), so it is converted to frequency for the FRO frequency axis using
# the rest frequency in the header.
#
# Author: Francesco Di Giovanni / HALO project
# Date: 2026-07-31

import os
import sys
import subprocess
import numpy as np
from datetime import datetime, timezone
from astropy.io import fits
import warnings
from astropy.wcs import WCS, FITSFixedWarning
from astropy.coordinates import SkyCoord, Galactic, FK5
import astropy.units as u

sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
import fro_format_v17 as fro

C_MPS = 299792458.0  # m/s

# JCMT (James Clerk Maxwell Telescope), Maunakea, Hawaii
JCMT_LAT = 19.8228
JCMT_LON = -155.4770
JCMT_ALT = 4092.0


def get_fits_path():
    """Return the path to the COHRS FITS cube: from argv[1] if given,
    otherwise ask via zenity file dialog. Never hardcoded to a specific tile."""
    if len(sys.argv) > 1:
        return sys.argv[1]
    surveys_dir = os.path.expanduser('~/FRO/Parsers/JCMT/Surveys/')
    result = subprocess.run(
        ['zenity', '--file-selection',
         '--title=Select COHRS FITS cube',
         '--file-filter=FITS files (*.fit *.fits)|*.fit *.fits',
         f'--filename={surveys_dir}'],
        capture_output=True, text=True)
    path = result.stdout.strip()
    if not path:
        print('No file selected, exiting.')
        sys.exit(1)
    return path


def load_cube(filepath):
    with fits.open(filepath) as hdul:
        data = np.asarray(hdul[0].data, dtype=np.float32)
        header = hdul[0].header.copy()
    return data, header


def build_velocity_axis_mps(header):
    """Reconstruct the LSRK radio velocity axis (m/s) from the FITS header.
    Nothing hardcoded: CRPIX3/CRVAL3/CDELT3/CUNIT3 read dynamically.
    COHRS uses CUNIT3=km/s, so the axis is scaled to m/s here."""
    n = header['NAXIS3']
    crpix3 = header['CRPIX3']
    crval3 = header['CRVAL3']
    cdelt3 = header['CDELT3']
    pix = np.arange(1, n + 1)
    v = crval3 + (pix - crpix3) * cdelt3
    cunit = str(header.get('CUNIT3', 'km/s')).strip().lower()
    if cunit in ('km/s', 'kms-1', 'km s-1'):
        v = v * 1000.0          # -> m/s
    # if already m/s, leave as is
    return v


def rest_frequency_hz(header):
    for k in ('RESTFRQ', 'RESTFREQ'):
        if k in header:
            return float(header[k])
    raise KeyError('No RESTFRQ/RESTFREQ in header; cannot build frequency axis.')


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


def convert_point(data, header, wcs2d, glon_req, glat_req, output_dir,
                  tile_name, timestamp):
    restfreq_hz = rest_frequency_hz(header)

    spectrum, ix, iy, glon_actual, glat_actual = extract_spectrum(
        data, header, wcs2d, glon_req, glat_req)

    v_mps = build_velocity_axis_mps(header)
    freq_hz = velocity_to_freq(v_mps, restfreq_hz)

    sort_idx = np.argsort(freq_hz)
    freq_hz_sorted = freq_hz[sort_idx]
    spectrum_sorted = spectrum[sort_idx]

    gc = SkyCoord(l=glon_actual * u.deg, b=glat_actual * u.deg, frame=Galactic())
    eq = gc.transform_to(FK5(equinox='J2000'))
    ra_deg = eq.ra.deg
    dec_deg = eq.dec.deg

    molecule = str(header.get('MOLECULE', '12CO(3-2)')).strip()

    out_name = f'COHRS_G{glon_actual:07.3f}{glat_actual:+07.3f}_FRO_{timestamp}.fro'
    out_path = os.path.join(output_dir, out_name)

    notes = (
        f'COHRS 12CO(3-2) High-Resolution Survey of the Galactic Plane '
        f'(Park et al. 2023, ApJS 264, 16; JCMT/HARP). '
        f'Source tile: {tile_name}. '
        f'Requested coordinates: GLon={glon_req:.4f}, GLat={glat_req:.4f}. '
        f'Actual extracted pixel-center coordinates: GLon={glon_actual:.4f}, '
        f'GLat={glat_actual:.4f} (pixel {ix},{iy}; grid step '
        f'{abs(header["CDELT1"]):.5f} deg). '
        f'Frequency axis reconstructed from the cube velocity axis '
        f'(CTYPE3={header.get("CTYPE3","VRAD")}, CUNIT3={header.get("CUNIT3","km/s")}, '
        f'LSRK) via f = f0*(1-v/c), f0 = RESTFRQ = {restfreq_hz:.2f} Hz '
        f'(from FITS header). Line: {molecule}. '
        f'Observation epoch: not defined for this product. COHRS is a '
        f'basket-weave OTF mosaic combining JCMT observations from several '
        f'semesters (project codes M10AU20, M11AU12, M11AD02, M12AD02, '
        f'M12AU42, M13AU41, M13BN02, M14AU09, M17BL004), rebinned onto a '
        f'common 0.635 km/s grid; Time/utc is therefore set to the explicit '
        f'UNKNOWN_EPOCH sentinel rather than to any date.'
    )

    fro.create_fro_v17_structure(
        filename=out_path,
        frequency_axis_hz=freq_hz_sorted,
        facility='JCMT 15m',
        instrument='HARP',
        observation_id=f'{tile_name}_px{ix}_{iy}',
        survey='COHRS',
        dataset=f'COHRS GLon{glon_actual:.3f} GLat{glat_actual:.3f}',
        provider='COHRS team (Park et al. 2023); JCMT / CADC',
        import_date=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        observatory_name='JCMT 15m',
        latitude_deg=JCMT_LAT,
        longitude_deg=JCMT_LON,
        elevation_m=JCMT_ALT,
        receiver='HARP (325-375 GHz)',
        center_frequency_hz=restfreq_hz,
        notes=notes,
        pointing_source='measured',
        temperature_calibration_note=(
            'Antenna temperature T_A* (K) on the COHRS calibrated scale '
            '(BUNIT=K in source FITS). Main-beam conversion factor ~0.61 '
            '(from survey README) to obtain T_mb; not applied here. '
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
        polarization='unknown',
        flag_track=-1,
        flag_cal=0,
    )

    return out_path


def main():
    print('=' * 60)
    print('  COHRS to FRO converter - HALO project')
    print('  Interactive point extraction from position-position-velocity cube')
    print('=' * 60)

    fits_path = os.path.expanduser(get_fits_path())
    if not os.path.isfile(fits_path):
        print(f'File not found: {fits_path}')
        sys.exit(1)

    tile_name = os.path.splitext(os.path.basename(fits_path))[0]

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

    print(f'Cube shape (vel, glat, glon): {data.shape}')
    print(f'RESTFRQ: {rest_frequency_hz(header):.2f} Hz')
    print(f'Grid step: {abs(header["CDELT1"]):.5f} deg')

    # Epoch directory derived from cube name + current timestamp.
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    surveys_root = os.path.join(os.path.dirname(fits_path), tile_name)
    epoch_dir = os.path.join(surveys_root, f'{tile_name}_FRO_{timestamp}')
    fro_dir   = os.path.join(epoch_dir, 'FRO')
    plots_dir = os.path.join(epoch_dir, 'Plots')
    raw_dir   = os.path.join(epoch_dir, 'Raw')
    os.makedirs(fro_dir,   exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)
    os.makedirs(raw_dir,   exist_ok=True)

    # Write source.txt in Raw/ pointing to the original cube — never copy it.
    source_txt = os.path.join(raw_dir, 'source.txt')
    with open(source_txt, 'w') as f:
        f.write(f'Source cube: {fits_path}\n')
    print(f'\nEpoch directory: {epoch_dir}')
    print(f'source.txt written in Raw/ -> {fits_path}')
    print(f'Output directory: {fro_dir}\n')

    while True:
        try:
            glon_req = float(input('GLon (deg): '))
            glat_req = float(input('GLat (deg): '))
        except ValueError:
            print('Invalid input, try again.\n')
            continue

        try:
            out_path = convert_point(
                data, header, wcs2d, glon_req, glat_req, fro_dir,
                tile_name, timestamp)
            print(f'  OK -> {os.path.basename(out_path)}')
        except Exception as e:
            print(f'  ERROR: {e}')

        again = input('\nExtract another point? [y/N]: ').strip().lower()
        if again != 'y':
            break

    print('\nDone.')

    # Hand over to the plotting suite, on the same cube, if wanted.
    plotter = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'view_cohrs.py')
    if os.path.isfile(plotter):
        if input('\nGenerate plots now? [y/N]: ').strip().lower() == 'y':
            subprocess.call([sys.executable, plotter, fits_path, epoch_dir])


if __name__ == '__main__':
    main()