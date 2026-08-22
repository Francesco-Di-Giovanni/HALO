# atca_to_fro.py
# Parser: ATCA (CASA-calibrated spectrum) → FRO HDF5 format v17
# Architecture: Option B — import calibrated product (spectrum .txt) from CASA.
# Raw interferometric visibilities are NOT stored; only the calibrated
# averaged spectrum is imported (consistent with the Option B decision
# documented in ATCA_Data_Structure_EN/IT_v1.3.docx).
# Reference dataset: project C1025, NGC 253 HI, 2002-02-09, config 1.5A.
# FRO/HALO Project — Francesco Di Giovanni (Claude AI assisted), Bolzano, Italy

import os
import sys
import numpy as np
from datetime import datetime, timezone
from astropy.coordinates import SkyCoord
import astropy.units as u

# CASA python must be on path — run with casa-python
try:
    from casatools import table as casatable
    CASA_AVAILABLE = True
except ImportError:
    CASA_AVAILABLE = False
    print('WARNING: casatools not available — metadata will use fallback values.')

sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
from fro_format_v17 import create_fro_v17_structure, append_raw_spectrum, read_fro_summary

# ── Paths ─────────────────────────────────────────────────────────────────────
_base       = os.path.expanduser('~/FRO/Parsers/ATCA/')
MS_PATH     = os.path.join(_base, 'Surveys/ngc253.ms')
SPEC_TXT    = os.path.join(_base, 'Surveys/NGC253/ngc253_spectrum.txt')
OUTPUT_DIR  = os.path.join(_base, 'Surveys/NGC253/')
OUTPUT_NAME = 'NGC253_ATCA_C1025_HI.fro'
OUTPUT_PATH = os.path.join(OUTPUT_DIR, OUTPUT_NAME)
# ─────────────────────────────────────────────────────────────────────────────


def read_ms_metadata(ms_path):
    meta = {}
    if not CASA_AVAILABLE:
        meta['ra_deg']             = 11.8974
        meta['dec_deg']            = -25.2939
        meta['center_freq_hz']     = 1419.0e6
        meta['bandwidth_hz']       = 8.0e6
        meta['integration_time_s'] = 29.2
        meta['obs_date']           = '2002-02-09T23:42:10'
        return meta

    tb = casatable()

    tb.open(os.path.join(ms_path, 'FIELD'))
    phase_dir = tb.getcell('PHASE_DIR', 2)
    ra_rad  = float(phase_dir[0]) if phase_dir.ndim == 1 else float(phase_dir[0, 0])
    dec_rad = float(phase_dir[1]) if phase_dir.ndim == 1 else float(phase_dir[1, 0])
    meta['ra_deg']  = np.degrees(ra_rad) % 360.0
    meta['dec_deg'] = np.degrees(dec_rad)
    tb.close()

    tb.open(os.path.join(ms_path, 'SPECTRAL_WINDOW'))
    meta['center_freq_hz'] = float(tb.getcell('REF_FREQUENCY', 0))
    meta['bandwidth_hz']   = float(tb.getcell('TOTAL_BANDWIDTH', 0))
    tb.close()

    tb.open(ms_path)
    field_col = tb.getcol('FIELD_ID')
    time_col  = tb.getcol('TIME')
    intv_col  = tb.getcol('INTERVAL')
    ngc_mask  = field_col == 2
    if np.any(ngc_mask):
        t_start = float(np.min(time_col[ngc_mask]))
        meta['obs_date']           = _mjd_s_to_utc(t_start).isoformat()
        meta['integration_time_s'] = float(np.median(intv_col[ngc_mask]))
    else:
        meta['obs_date']           = '2002-02-09T23:42:10'
        meta['integration_time_s'] = 29.2
    tb.close()

    return meta


def _mjd_s_to_utc(mjd_s):
    MJD_UNIX_OFFSET = 3506716800.0
    unix_s = mjd_s - MJD_UNIX_OFFSET
    return datetime.fromtimestamp(unix_s, tz=timezone.utc)


def read_spectrum(spec_txt):
    data = np.loadtxt(spec_txt)
    freq_hz   = data[:, 0] * 1e6
    amplitude = data[:, 1]
    return freq_hz, amplitude


def compute_galactic(ra_deg, dec_deg):
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    return coord.galactic.l.deg, coord.galactic.b.deg


def build_fro(ms_path, spec_txt, output_path):
    print('Reading MS metadata from:', ms_path)
    meta = read_ms_metadata(ms_path)

    print('Reading calibrated spectrum from:', spec_txt)
    freq_hz, amplitude = read_spectrum(spec_txt)

    n_chan = len(freq_hz)
    print(f'Spectrum: {n_chan} channels, {freq_hz[0]/1e6:.3f}–{freq_hz[-1]/1e6:.3f} MHz')
    print(f'Amplitude range: {amplitude.min():.4f}–{amplitude.max():.4f} Jy')

    ra_deg  = meta['ra_deg']
    dec_deg = meta['dec_deg']
    glon, glat = compute_galactic(ra_deg, dec_deg)

    obs_utc_str = meta['obs_date']
    try:
        obs_utc = datetime.fromisoformat(obs_utc_str)
        if obs_utc.tzinfo is None:
            obs_utc = obs_utc.replace(tzinfo=timezone.utc)
    except ValueError:
        obs_utc = datetime(2002, 2, 9, 23, 42, 10, tzinfo=timezone.utc)

    notes = (
        'ATCA project C1025 (PI: R. Boomsma). Target: NGC 253. '
        'Observation date: 2002-02-09. Array configuration: 1.5A. '
        'Architecture: Option B — calibrated spectrum imported from CASA 6.7.5-18. '
        'Calibration: bandpass from 1934-638 (scans 0-1), '
        'gain/phase from 0023-263 (scan 2), applied to NGC 253 (scan 3). '
        'Spectrum extracted from CORRECTED_DATA column, SpW 0 (IF0, HI band), '
        'averaged over all 390 visibilities (15 baselines x 26 integrations). '
        'Amplitude spectrum (not flux density): raw interferometric amplitude '
        'after calibration, not deconvolved or imaged. '
        'Continuum NOT subtracted. '
        'Pointing coordinates from FIELD table, field id 2 (NGC 253). '
        'Galactic coordinates computed from RA/Dec via Astropy. '
        'Raw RPFITS file: 2002-02-09_2302.C1025 (ATOA archive, OPAL). '
        'CASA Measurement Set: ngc253.ms. '
        'Calibration tables: bandpass.cal, gain.cal.'
    )

    create_fro_v17_structure(
        output_path,
        frequency_axis_hz  = freq_hz,
        facility           = 'ATCA',
        instrument         = 'CA01-CA06 (6x22m, config 1.5A)',
        observation_id     = 'C1025',
        survey             = 'ATCA archival C1025',
        dataset            = 'NGC253_HI_2002-02-09',
        provider           = 'ATNF / ATOA',
        import_date        = datetime.now(timezone.utc).isoformat(),
        observatory_name   = 'Australia Telescope Compact Array (ATCA)',
        latitude_deg       = -30.3128,
        longitude_deg      = 149.5501,
        elevation_m        = 237.0,
        receiver           = 'L-band (1.4 GHz)',
        antenna            = '6 x 22m antennas',
        lna                = 'unknown',
        center_frequency_hz = meta.get('center_freq_hz', 1419.0e6),
        bandwidth_hz       = meta.get('bandwidth_hz', 8.0e6),
        integration_time_s = meta.get('integration_time_s', 29.2),
        pointing_source    = 'measured',
        target             = 'NGC 253',
        notes              = notes,
        temperature_calibration_note = (
            'Amplitude in Jy (interferometric visibility amplitude after '
            'bandpass and gain/phase calibration via CASA). Not a single-dish '
            'antenna temperature. Continuum not subtracted.'
        ),
    )
    print('Created FRO file:', output_path)

    append_raw_spectrum(
        output_path,
        raw_spectrum       = amplitude.astype(np.float32),
        utc_time           = obs_utc,
        integration_time_s = meta.get('integration_time_s', 29.2),
        ra_deg             = ra_deg,
        dec_deg            = dec_deg,
        glon_deg           = glon,
        glat_deg           = glat,
        az_deg             = np.nan,
        elev_deg           = np.nan,
        polarization       = 'XX+YY',
        flag_track         = 1,
        flag_cal           = 0,
    )
    print('Spectrum written (1 row, XX+YY averaged).')
    print()
    read_fro_summary(output_path)


if __name__ == '__main__':
    build_fro(MS_PATH, SPEC_TXT, OUTPUT_PATH)
