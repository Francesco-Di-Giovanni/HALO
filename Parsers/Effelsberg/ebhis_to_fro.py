# ebhis_to_fro.py
# Converts EBHIS survey .txt profiles to FRO HDF5 format (v17)
# EBHIS = Effelsberg-Bonn HI Survey (100m Effelsberg telescope)
# Data source: https://www.astro.uni-bonn.de/hisurvey/AllSky_profiles/
#
# Frequency axis: converted from v_lsr [km/s] to frequency [Hz]
# using f = f0 * (1 - v/c), f0 = 1420405750 Hz (HI rest frequency)
# Conversion documented in Notes field of each .fro file.
#
# Author: Francesco Di Giovanni / HALO project
# Date: 2026-07-22

import os
import sys
import numpy as np
from datetime import datetime, timezone
from astropy.coordinates import SkyCoord, Galactic, FK5
import astropy.units as u

# Add FRO_System to path for fro_format_v17
sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
import fro_format_v17 as fro

# Constants
HI_REST_FREQ_HZ = 1420405750.0  # Hz
C_KMS = 299792.458              # km/s

# ── Grid settings — asked interactively at runtime ─────────────────────────────
GLON_MIN = float(input('GLon min (deg): '))
GLON_MAX = float(input('GLon max (deg): '))
GLAT_MIN = float(input('GLat min (deg): '))
GLAT_MAX = float(input('GLat max (deg): '))

# Auto-build directory names from grid coordinates
GRID_TAG   = f'G{GLON_MIN:03.0f}-{GLON_MAX:03.0f}_B{GLAT_MIN:+03.0f}{GLAT_MAX:+03.0f}'
SURVEY_DIR = os.path.expanduser(f'~/FRO/Parsers/Effelsberg/Surveys/EBHIS_Survey_{GRID_TAG}')
EBHIS_DIR  = os.path.join(SURVEY_DIR, f'EBHIS_Profiles_{GRID_TAG}')
OUTPUT_DIR = os.path.join(SURVEY_DIR, f'EBHIS_FRO_{GRID_TAG}')
# ───────────────────────────────────────────────────────────────────────────────

# Observatory: Effelsberg 100m
OBS_LAT = 50.5248    # deg N
OBS_LON = 6.8836     # deg E
OBS_ALT = 369.0      # m


def parse_ebhis_file(filepath):
    """Parse a single EBHIS .txt profile file.
    Returns dict with metadata and spectral data."""
    glon = None
    glat = None
    n_points = None
    nh = None
    velocities = []
    tb_values = []

    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('%%') and 'datapoints' not in line:
                # First line: %% glon glat 0 beam
                parts = line.lstrip('%').split()
                glon = float(parts[0])
                glat = float(parts[1])
            elif 'datapoints' in line:
                # %%EBHIS   621  datapoints: ...
                parts = line.split()
                n_points = int(parts[1])
            elif line.startswith('%') and 'NH  within specified beam' in line:
                parts = line.split()
                nh = float(parts[-1])
            elif not line.startswith('%') and line.strip():
                parts = line.split()
                if len(parts) == 2:
                    velocities.append(float(parts[0]))
                    tb_values.append(float(parts[1]))

    velocities = np.array(velocities, dtype=np.float32)
    tb_values = np.array(tb_values, dtype=np.float32)

    # Convert v_lsr [km/s] -> frequency [Hz]
    freq_hz = HI_REST_FREQ_HZ * (1.0 - velocities / C_KMS)

    # Sort by increasing frequency
    sort_idx = np.argsort(freq_hz)
    freq_hz = freq_hz[sort_idx]
    tb_values = tb_values[sort_idx]

    # Compute RA/Dec from Galactic coordinates
    gc = SkyCoord(l=glon * u.deg, b=glat * u.deg, frame=Galactic())
    eq = gc.transform_to(FK5(equinox='J2000'))
    ra_deg = eq.ra.deg
    dec_deg = eq.dec.deg

    return {
        'glon': glon,
        'glat': glat,
        'ra_deg': ra_deg,
        'dec_deg': dec_deg,
        'nh': nh,
        'n_points': n_points,
        'freq_hz': freq_hz,
        'tb_k': tb_values,
    }



def get_files_in_range(glon_min, glon_max, glat_min, glat_max):
    """Return sorted list of EBHIS .txt files within the specified range."""
    files = []
    for fname in sorted(os.listdir(EBHIS_DIR)):
        if not fname.endswith('.txt'):
            continue
        try:
            core = fname.replace('EBHIS_G', '').replace('.txt', '')
            if '+' in core:
                lon_str, lat_str = core.split('+')
                lat_val = float(lat_str)
            else:
                parts = core.split('-')
                lon_str = parts[0]
                lat_val = -float(parts[1])
            lon_val = float(lon_str)
            if glon_min <= lon_val <= glon_max and glat_min <= lat_val <= glat_max:
                files.append(fname)
        except Exception:
            continue
    return files


def convert_file(fname, obs_date):
    """Convert a single EBHIS .txt file to .fro v17."""
    filepath = os.path.join(EBHIS_DIR, fname)
    data = parse_ebhis_file(filepath)

    glon = data['glon']
    glat = data['glat']
    freq_hz = data['freq_hz']
    tb_k = data['tb_k']

    # Output filename
    out_name = f"EBHIS_G{glon:06.2f}{glat:+06.2f}_FRO.fro"
    out_path = os.path.join(OUTPUT_DIR, out_name)

    notes = (
        f"EBHIS survey profile. Source file: {fname}. "
        f"Frequency axis converted from v_lsr [km/s] using f = f0*(1-v/c), "
        f"f0 = {HI_REST_FREQ_HZ:.0f} Hz (HI rest frequency). "
        f"T_B [K] = brightness temperature (beam-averaged). "
        f"NH = {data['nh']:.3e} cm^-2 (from file header). "
        f"N_datapoints = {data['n_points']}."
    )

    fro.create_fro_v17_structure(
        filename=out_path,
        frequency_axis_hz=freq_hz,
        facility='Effelsberg',
        instrument='Effelsberg 100m',
        observation_id=fname.replace('.txt', ''),
        survey='EBHIS',
        dataset=f'EBHIS GLon{glon:.2f} GLat{glat:.2f}',
        provider='MPIfR Bonn / Argelander Institut',
        import_date=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
        observatory_name='Effelsberg',
        latitude_deg=OBS_LAT,
        longitude_deg=OBS_LON,
        elevation_m=OBS_ALT,
        receiver='L-band',
        center_frequency_hz=HI_REST_FREQ_HZ,
        notes=notes,
        pointing_source='computed',
    )

    fro.append_raw_spectrum(
        filename=out_path,
        raw_spectrum=tb_k,
        utc_time=obs_date,
        az_deg=float('nan'),
        elev_deg=float('nan'),
        ra_deg=data['ra_deg'],
        dec_deg=data['dec_deg'],
        glon_deg=glon,
        glat_deg=glat,
        integration_time_s=0.0,
        polarization='I',
        flag_track=-1,
        flag_cal=0,
    )

    return out_path


def main():
    print('=' * 60)
    print('  EBHIS to FRO converter — HALO project')
    print('=' * 60)
    print(f'\nEBHIS source directory : {EBHIS_DIR}')
    print(f'Output directory       : {OUTPUT_DIR}')
    print()
    files = get_files_in_range(GLON_MIN, GLON_MAX, GLAT_MIN, GLAT_MAX)

    if not files:
        print('\nNo EBHIS files found in the specified range.')
        print(f'Available range: GLon {GLON_MIN}-{GLON_MAX}, GLat {GLAT_MIN}/+{GLAT_MAX}')
        sys.exit(0)

    print(f'\nFound {len(files)} files in the specified range.')
    print('First:', files[0])
    print('Last: ', files[-1])

    confirm = input('\nProceed with conversion? [y/N]: ').strip().lower()
    if confirm != 'y':
        print('Aborted.')
        sys.exit(0)

    # Use a fixed reference date for EBHIS (survey epoch ~2008-2013)
    obs_date = datetime(2010, 6, 15, 12, 0, 0, tzinfo=timezone.utc)

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print()
    errors = 0
    for i, fname in enumerate(files, 1):
        try:
            out_path = convert_file(fname, obs_date)
            print(f'  [{i:3d}/{len(files)}] OK  {os.path.basename(out_path)}')
        except Exception as e:
            print(f'  [{i:3d}/{len(files)}] ERR {fname}: {e}')
            errors += 1

    print()
    print(f'Done. {len(files) - errors} converted, {errors} errors.')
    print(f'Output: {OUTPUT_DIR}')


if __name__ == '__main__':
    main()