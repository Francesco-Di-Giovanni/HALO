# parkes_gass_to_fro.py
# Parser: Parkes GASS SDFITS -> FRO HDF5 format v17
# Parkes 64m multibeam receiver (7 beams, 2 IFs, 2 pols, 2048 channels)
# Reference file: Parkes_GASS.fits (NASA FITS registry sample)
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy

import os
import sys
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, Galactic, FK5, AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u
from datetime import datetime, timezone

sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
from fro_format_v17 import create_fro_v17_structure, append_raw_spectrum, read_fro_summary

# Parkes observatory coordinates
PARKES_LAT = -32.9984  # deg
PARKES_LON = 148.2635  # deg
PARKES_ALT = 415.0     # m
PARKES_LOCATION = EarthLocation(
    lat=PARKES_LAT * u.deg,
    lon=PARKES_LON * u.deg,
    height=PARKES_ALT * u.m
)

# Polarization map: STOKES CRVAL2 codes
STOKES_MAP = {1: 'I', 2: 'Q', 3: 'U', 4: 'V',
              -1: 'RR', -2: 'LL', -3: 'RL', -4: 'LR',
              -5: 'XX', -6: 'YY', -7: 'XY', -8: 'YX'}

INPUT_FILE = os.path.expanduser('~/FRO/Parsers/Parkes/Surveys/Parkes_Survey_GASS/Parkes_Profiles_GASS/Parkes_GASS.fits')
OUTPUT_DIR = os.path.expanduser('~/FRO/Parsers/Parkes/Surveys/Parkes_Survey_GASS/Parkes_FRO_GASS/')


def get_freq_axis(row):
    """Build frequency axis in Hz from WCS keywords."""
    crval = float(row['CRVAL1'])
    cdelt = float(row['CDELT1'])
    crpix = float(row['CRPIX1'])
    n = 2048
    channels = np.arange(1, n + 1, dtype=np.float64)
    return crval + (channels - crpix) * cdelt


def get_radec(row):
    """Extract RA/Dec from CRVAL3/CRVAL4."""
    ra = float(row['CRVAL3'])
    dec = float(row['CRVAL4'])
    return ra, dec


def get_galactic(ra_deg, dec_deg):
    """Convert RA/Dec to Galactic coordinates."""
    eq = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame=FK5(equinox='J2000'))
    gal = eq.galactic
    return float(gal.l.deg), float(gal.b.deg)


def get_altaz(ra_deg, dec_deg, date_obs, time_s):
    """Compute Az/El from RA/Dec and time."""
    try:
        t = Time(f'{date_obs}T00:00:00', format='isot', scale='utc') + time_s * u.s
        coord = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame=FK5(equinox='J2000'))
        altaz_frame = AltAz(obstime=t, location=PARKES_LOCATION)
        altaz = coord.transform_to(altaz_frame)
        return float(altaz.az.deg), float(altaz.alt.deg)
    except Exception:
        return float('nan'), float('nan')


def get_utc(date_obs, time_s):
    """Build UTC datetime from DATE-OBS and TIME."""
    try:
        t = Time(f'{date_obs}T00:00:00', format='isot', scale='utc') + float(time_s) * u.s
        return t.to_datetime(timezone=timezone.utc)
    except Exception:
        return datetime(2004, 9, 15, 0, 0, 0, tzinfo=timezone.utc)


def parse_parkes_gass(input_file, output_dir):
    print(f'Opening: {input_file}')
    h = fits.open(input_file)
    t = h[1]
    data_table = t.data
    n_rows = len(data_table)
    print(f'Rows: {n_rows}')

    # Group by BEAM — one .fro per beam
    beams = np.unique(data_table['BEAM'])
    print(f'Beams: {beams}')

    for beam in beams:
        mask = data_table['BEAM'] == beam
        rows = data_table[mask]
        n = np.sum(mask)
        print(f'\nBeam {beam}: {n} rows')

        # Frequency axis from first row
        freq_hz = get_freq_axis(rows[0])

        # Output filename
        out_path = os.path.join(output_dir, f'Parkes_GASS_beam{beam:02d}.fro')

        notes = (
            f'Parkes GASS SDFITS sample. Beam {beam}. '
            f'Frequency-switched data (OBSMODE={str(rows[0]["OBSMODE"]).strip()}). '
            f'CRVAL1 center freq offset from HI rest frequency (frequency switching). '
            f'TSYS vector: 2 values per row (one per polarization). '
            f'TAMBIENT/PRESSURE/HUMIDITY not recorded in this sample (all zero). '
            f'Az/El computed from RA/Dec via Astropy.'
        )

        create_fro_v17_structure(
            filename=out_path,
            frequency_axis_hz=freq_hz,
            facility='Parkes',
            instrument='Parkes 64m multibeam',
            observation_id=f'GASS_beam{beam:02d}',
            survey='GASS',
            dataset=f'Parkes GASS beam {beam}',
            provider='CSIRO ATNF',
            import_date=datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ'),
            observatory_name='Parkes Observatory',
            latitude_deg=PARKES_LAT,
            longitude_deg=PARKES_LON,
            elevation_m=PARKES_ALT,
            receiver=f'Multibeam 21cm beam {beam}',
            center_frequency_hz=1420405750.0,
            notes=notes,
            pointing_source='computed',
        )

        for i, row in enumerate(rows):
            ra, dec = get_radec(row)
            glon, glat = get_galactic(ra, dec)
            az, el = get_altaz(ra, dec, row['DATE-OBS'], row['TIME'])
            utc = get_utc(row['DATE-OBS'], row['TIME'])
            integration = float(row['EXPOSURE'])

            # DATA shape: (1, 1, 2, 2048) — pols x chans
            raw_data = row['DATA'].squeeze()  # -> (2, 2048)
            tsys = row['TSYS']  # [tsys_pol0, tsys_pol1]

            # Polarization labels from STOKES header
            pol_labels = ['XX', 'YY']  # Parkes multibeam uses linear feeds

            for p in range(2):
                spectrum = raw_data[p, :]  # 2048 channels

                append_raw_spectrum(
                    filename=out_path,
                    raw_spectrum=spectrum.astype(np.float32),
                    utc_time=utc,
                    az_deg=az,
                    elev_deg=el,
                    ra_deg=ra,
                    dec_deg=dec,
                    glon_deg=glon,
                    glat_deg=glat,
                    integration_time_s=integration,
                    polarization=pol_labels[p],
                    flag_track=0,
                    flag_cal=0,
                    temperature_c=float(row['TAMBIENT']) if row['TAMBIENT'] != 0 else float('nan'),
                    pressure_mbar=float(row['PRESSURE']) if row['PRESSURE'] != 0 else float('nan'),
                    humidity_percent=float(row['HUMIDITY']) * 100.0 if row['HUMIDITY'] != 0 else float('nan'),
                )

        print(f'  -> {out_path}')
        read_fro_summary(out_path)

    h.close()
    print('\nDone.')


if __name__ == '__main__':
    if len(sys.argv) > 1:
        parse_parkes_gass(sys.argv[1], OUTPUT_DIR)
    else:
        parse_parkes_gass(INPUT_FILE, OUTPUT_DIR)