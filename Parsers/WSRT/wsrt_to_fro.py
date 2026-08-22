#!/usr/bin/env python3
"""
wsrt_to_fro.py — WSRT / Apertif to HALO v17 parser
Supports two input modes:
  - MeasurementSet (.ms directory): extracts autocorrelations via casatools
  - FITS cube / UVFITS (.fits): extracts spectra via astropy

For MeasurementSet input, autocorrelations (ANTENNA1 == ANTENNA2) are
extracted as single-dish total-power spectra per antenna and integration.
For FITS cube input, an integrated spectrum is computed over all spatial axes.
UVFITS input is detected via SIMPLE/GROUPS header and handled accordingly.

Usage:
    python3 wsrt_to_fro.py [input_path]

Requires:
  - astropy, numpy, h5py (standard HALO environment)
  - casatools (conda env 'casa_env') for .ms input mode
  - fro_format_v17.py in ~/FRO/FRO_System/

Output:
  ~/FRO/Parsers/WSRT/Surveys/<dataset>_<YYYYMMDD_HHMMSS>/FRO/<name>.fro

Project: HALO (Hydrogen Atomic Line Observatory)
Author: Francesco Di Giovanni (Bolzano, Italy)
Created: 2026-08-21 (Parte 67)
"""

import os
import sys
import subprocess
import numpy as np
from datetime import datetime, timezone
from pathlib import Path

# ── HALO system path ─────────────────────────────────────────────────────────
FRO_SYSTEM = Path.home() / 'FRO' / 'FRO_System'
sys.path.insert(0, str(FRO_SYSTEM))
import fro_format_v17 as fro

# ── Output base directory ────────────────────────────────────────────────────
SURVEYS_DIR = Path.home() / 'FRO' / 'Parsers' / 'WSRT' / 'Surveys'

# ── WSRT observatory coordinates (Westerbork, Netherlands) ──────────────────
WSRT_LAT_DEG  =  52.9142
WSRT_LON_DEG  =   6.6033
WSRT_ELEV_M   =  16.0

# ── Sentinel for unknown epoch ───────────────────────────────────────────────
UNKNOWN_EPOCH = fro.UNKNOWN_EPOCH


def select_input():
    """Prompt user to select input file or directory via zenity."""
    try:
        result = subprocess.run(
            ['zenity', '--file-selection',
             '--title=Seleziona MeasurementSet (.ms) o FITS/UVFITS',
             '--filename=' + str(Path.home() / 'Downloads') + '/',
             '--file-filter=FITS files | *.fits *.fit'],
            capture_output=True, text=True
        )
        path = result.stdout.strip()
        if not path:
            print("Nessun file selezionato. Uscita.")
            sys.exit(0)
        return Path(path)
    except FileNotFoundError:
        path = input("Percorso MeasurementSet o FITS: ").strip()
        if not path:
            sys.exit(0)
        return Path(path)


def detect_mode(input_path):
    """Return 'ms', 'uvfits', or 'fits' based on input."""
    p = str(input_path).lower()
    if input_path.is_dir() and (p.endswith('.ms') or (input_path / 'table.info').exists()):
        return 'ms'
    if input_path.is_dir():
        return 'ms'
    if input_path.is_file() and (p.endswith('.fits') or p.endswith('.fit')):
        # Distinguish UVFITS (GROUPS=T) from image FITS cube
        from astropy.io import fits as astrofits
        with astrofits.open(str(input_path)) as hdul:
            hdr = hdul[0].header
            if hdr.get('GROUPS', False):
                return 'uvfits'
        return 'fits'
    print(f"Tipo input non riconosciuto: {input_path}")
    sys.exit(1)


# ════════════════════════════════════════════════════════════════════════════
# Mode A: MeasurementSet via casatools
# ════════════════════════════════════════════════════════════════════════════

def parse_ms(ms_path):
    """Extract autocorrelation spectra from a WSRT MeasurementSet."""
    print(f"Modalità: MeasurementSet")
    print(f"Input: {ms_path}")

    try:
        from casatools import ms as mstool, table as tbtool
    except ImportError:
        print("\nErrore: casatools non trovato.")
        print("Attivare l'ambiente conda: conda activate casa_env")
        sys.exit(1)

    tb = tbtool()

    tb.open(str(ms_path) + '/ANTENNA')
    antenna_names = list(tb.getcol('NAME'))
    tb.close()

    tb.open(str(ms_path) + '/SPECTRAL_WINDOW')
    chan_freq = tb.getcol('CHAN_FREQ')
    ref_freq  = tb.getcol('REF_FREQUENCY')
    tb.close()

    tb.open(str(ms_path) + '/FIELD')
    phase_dir = tb.getcol('PHASE_DIR')
    tb.close()

    tb.open(str(ms_path))
    ant1     = tb.getcol('ANTENNA1')
    ant2     = tb.getcol('ANTENNA2')
    data     = tb.getcol('DATA')
    time_col = tb.getcol('TIME')
    spw_id   = tb.getcol('DATA_DESC_ID')
    tb.close()

    auto_mask = (ant1 == ant2)
    if not np.any(auto_mask):
        print("Nessuna autocorrelazione trovata.")
        sys.exit(1)

    data_auto  = data[:, :, auto_mask]
    time_auto  = time_col[auto_mask]
    spw_auto   = spw_id[auto_mask]

    n_auto    = data_auto.shape[2]
    spw_idx   = int(spw_auto[0])
    frequency_axis_hz   = chan_freq[:, spw_idx].astype(np.float64)
    center_hz = float(ref_freq[spw_idx])

    unique_times  = np.unique(time_auto)
    integration_s = float(np.median(np.diff(unique_times))) if len(unique_times) > 1 else float('nan')

    ra_deg  = float(np.degrees(phase_dir[0, 0, 0]))
    dec_deg = float(np.degrees(phase_dir[1, 0, 0]))

    from astropy.coordinates import SkyCoord
    import astropy.units as u
    coord    = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    gal      = coord.galactic
    glon_deg = float(gal.l.deg)
    glat_deg = float(gal.b.deg)

    dataset_name = ms_path.stem
    epoch_str    = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    epoch_dir    = SURVEYS_DIR / f"{dataset_name}_{epoch_str}"
    fro_dir      = epoch_dir / 'FRO'
    raw_dir      = epoch_dir / 'Raw'
    fro_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    with open(raw_dir / 'source.txt', 'w') as f:
        f.write(f"MeasurementSet: {ms_path}\n")
        f.write(f"Antennas: {', '.join(antenna_names)}\n")

    fro_path = fro_dir / f"{dataset_name}_{epoch_str}.fro"

    fro.create_fro_v17_structure(
        filename           = str(fro_path),
        facility           = 'WSRT',
        instrument         = 'WSRT / Apertif',
        survey             = dataset_name,
        dataset            = dataset_name,
        provider           = 'ASTRON',
        observation_id     = dataset_name,
        observatory_name   = 'WSRT',
        latitude_deg       = WSRT_LAT_DEG,
        longitude_deg      = WSRT_LON_DEG,
        elevation_m        = WSRT_ELEV_M,
        frequency_axis_hz  = frequency_axis_hz,
        center_frequency_hz= center_hz,
        receiver           = 'WSRT L-band',
        integration_time_s = n_auto,
        pointing_source    = 'measured',
        notes              = f"WSRT autocorrelations from {ms_path.name}. "
                             f"ANTENNA1==ANTENNA2 rows only. "
                             f"Parsed by wsrt_to_fro.py.",
        overwrite          = False
    )

    n_corr = data_auto.shape[0]
    if n_corr >= 4:
        stokes_i = 0.5 * (np.abs(data_auto[0, :, :]) + np.abs(data_auto[3, :, :]))
    else:
        stokes_i = np.abs(data_auto[0, :, :])

    for i in range(n_auto):
        from astropy.time import Time
        t      = Time(time_auto[i] / 86400.0, format='mjd', scale='utc')
        utc_dt = t.to_datetime(timezone.utc)

        fro.append_raw_spectrum(
            filename           = str(fro_path),
            raw_spectrum       = stokes_i[:, i].astype(np.float32),
            utc_time           = utc_dt,
            integration_time_s = integration_s,
            az_deg             = float('nan'),
            elev_deg           = float('nan'),
            ra_deg             = ra_deg,
            dec_deg            = dec_deg,
            glon_deg           = glon_deg,
            glat_deg           = glat_deg,
        )
        if (i + 1) % 100 == 0 or i == n_auto - 1:
            print(f"  Spettro {i+1}/{n_auto}", end='\r')

    print(f"\nFRO scritto: {fro_path}")
    return fro_path


# ════════════════════════════════════════════════════════════════════════════
# Mode B: UVFITS
# ════════════════════════════════════════════════════════════════════════════

def parse_uvfits(fits_path):
    """
    Extract autocorrelations from a WSRT UVFITS file.
    UVFITS stores visibilities as random groups; baseline 0 = autocorrelation.
    Baseline encoding: baseline = ant1 * 256 + ant2 (AIPS convention).
    Autocorrelations: ant1 == ant2, i.e. baseline % 257 == 0.
    """
    print(f"Modalità: UVFITS")
    print(f"Input: {fits_path}")

    from astropy.io import fits as astrofits
    import warnings
    from astropy.utils.exceptions import AstropyWarning

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        with astrofits.open(str(fits_path)) as hdul:
            hdr  = hdul[0].header
            data = hdul[0].data   # random groups

    # Frequency axis from header (FITS axis 4 = FREQ in UVFITS convention)
    nfreq  = int(hdr.get('NAXIS4', 1))
    crpix  = float(hdr.get('CRPIX4', 1.0))
    crval  = float(hdr.get('CRVAL4', 0.0))
    cdelt  = float(hdr.get('CDELT4', 1.0))
    frequency_axis_hz   = (crval + (np.arange(nfreq) - (crpix - 1)) * cdelt).astype(np.float64)
    center_hz = float(crval)

    # Time from DATE (JD) parameter
    date_obs = hdr.get('DATE-OBS', None)
    if date_obs:
        try:
            from astropy.time import Time
            utc_dt = Time(date_obs, format='isot', scale='utc').to_datetime(timezone.utc)
        except Exception:
            utc_dt = UNKNOWN_EPOCH
    else:
        utc_dt = UNKNOWN_EPOCH

    # Phase centre
    try:
        ra_deg  = float(hdr.get('CRVAL6', float('nan')))  # RA in degrees
        dec_deg = float(hdr.get('CRVAL7', float('nan')))  # Dec in degrees
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        coord    = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
        gal      = coord.galactic
        glon_deg = float(gal.l.deg)
        glat_deg = float(gal.b.deg)
    except Exception:
        ra_deg = dec_deg = glon_deg = glat_deg = float('nan')

    # Extract autocorrelations: BASELINE parameter, AIPS convention
    baselines = data.par('BASELINE')
    ant1 = (baselines / 256).astype(int)
    ant2 = (baselines % 256).astype(int)
    auto_mask = (ant1 == ant2)

    if not np.any(auto_mask):
        print("Nessuna autocorrelazione trovata nel UVFITS.")
        print("Procedendo con spettro integrato su tutte le visibilità.")
        auto_mask = np.ones(len(baselines), dtype=bool)

    # Visibility data shape: (nrows, 1, 1, nfreq, nstokes, 2) — real/imag
    vis = data.data[auto_mask]
    # Stokes I: mean of XX (0) and YY (1) real parts, weight > 0
    weights = vis[..., 2]   # weight is the 3rd element (real, imag, weight)
    real    = vis[..., 0]
    # Average over autocorrelation rows, stokes I approximation
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        raw_spectrum = np.nanmean(real[:, 0, 0, :, 0, ...], axis=0).astype(np.float32)
    if spectrum.ndim > 1:
        spectrum = spectrum.mean(axis=-1)

    dataset_name = fits_path.stem
    epoch_str    = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    epoch_dir    = SURVEYS_DIR / f"{dataset_name}_{epoch_str}"
    fro_dir      = epoch_dir / 'FRO'
    raw_dir      = epoch_dir / 'Raw'
    fro_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    with open(raw_dir / 'source.txt', 'w') as f:
        f.write(f"UVFITS: {fits_path}\n")

    fro_path = fro_dir / f"{dataset_name}_{epoch_str}.fro"

    fro.create_fro_v17_structure(
        filename           = str(fro_path),
        facility           = 'WSRT',
        instrument         = hdr.get('INSTRUME', 'WSRT'),
        survey             = hdr.get('OBJECT', dataset_name),
        dataset            = dataset_name,
        provider           = 'ASTRON',
        observation_id     = dataset_name,
        observatory_name   = 'WSRT',
        latitude_deg       = WSRT_LAT_DEG,
        longitude_deg      = WSRT_LON_DEG,
        elevation_m        = WSRT_ELEV_M,
        frequency_axis_hz  = frequency_axis_hz,
        center_frequency_hz= center_hz,
        receiver           = hdr.get('INSTRUME', 'WSRT'),
        pointing_source    = 'measured',
        notes              = f"WSRT UVFITS {fits_path.name}. "
                             f"Autocorrelation rows averaged to single spectrum. "
                             f"Parsed by wsrt_to_fro.py.",
        overwrite          = False
    )

    fro.append_raw_spectrum(
        filename           = str(fro_path),
        raw_spectrum       = raw_spectrum,
        utc_time           = utc_dt,
        integration_time_s = float('nan'),
        az_deg             = float('nan'),
        elev_deg           = float('nan'),
        ra_deg             = ra_deg,
        dec_deg            = dec_deg,
        glon_deg           = glon_deg,
        glat_deg           = glat_deg,
    )

    print(f"FRO scritto: {fro_path}")
    return fro_path


# ════════════════════════════════════════════════════════════════════════════
# Mode C: FITS cube (Apertif survey products)
# ════════════════════════════════════════════════════════════════════════════

def parse_fits(fits_path):
    """Extract integrated spectrum from a WSRT/Apertif FITS cube."""
    print(f"Modalità: FITS cube")
    print(f"Input: {fits_path}")

    from astropy.io import fits as astrofits
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    import warnings
    from astropy.utils.exceptions import AstropyWarning

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        with astrofits.open(str(fits_path)) as hdul:
            hdr  = hdul[0].header
            data = hdul[0].data

    n_axes    = data.ndim
    freq_axis = None
    for i in range(1, n_axes + 1):
        if 'FREQ' in hdr.get(f'CTYPE{i}', '').upper():
            freq_axis = i - 1
            break

    if freq_axis is None:
        print("Asse FREQ non trovato.")
        sys.exit(1)

    n_freq = data.shape[n_axes - 1 - freq_axis]
    crpix  = hdr.get(f'CRPIX{freq_axis+1}', 1.0)
    crval  = hdr.get(f'CRVAL{freq_axis+1}', 0.0)
    cdelt  = hdr.get(f'CDELT{freq_axis+1}', 1.0)
    frequency_axis_hz   = (crval + (np.arange(n_freq) - (crpix - 1)) * cdelt).astype(np.float64)
    center_hz = float(crval)

    date_obs = hdr.get('DATE-OBS', None)
    if date_obs:
        try:
            from astropy.time import Time
            utc_dt = Time(date_obs, format='isot', scale='utc').to_datetime(timezone.utc)
        except Exception:
            utc_dt = UNKNOWN_EPOCH
    else:
        utc_dt = UNKNOWN_EPOCH

    try:
        ra_deg   = float(hdr.get('CRVAL1', float('nan')))
        dec_deg  = float(hdr.get('CRVAL2', float('nan')))
        coord    = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
        gal      = coord.galactic
        glon_deg = float(gal.l.deg)
        glat_deg = float(gal.b.deg)
    except Exception:
        ra_deg = dec_deg = glon_deg = glat_deg = float('nan')

    freq_numpy_axis = n_axes - 1 - freq_axis
    axes_to_mean    = tuple(i for i in range(n_axes) if i != freq_numpy_axis)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        raw_spectrum = np.nanmean(data, axis=axes_to_mean).astype(np.float32)

    dataset_name = fits_path.stem
    epoch_str    = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    epoch_dir    = SURVEYS_DIR / f"{dataset_name}_{epoch_str}"
    fro_dir      = epoch_dir / 'FRO'
    raw_dir      = epoch_dir / 'Raw'
    fro_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    with open(raw_dir / 'source.txt', 'w') as f:
        f.write(f"FITS cube: {fits_path}\n")

    fro_path = fro_dir / f"{dataset_name}_{epoch_str}.fro"

    fro.create_fro_v17_structure(
        filename           = str(fro_path),
        facility           = 'WSRT',
        instrument         = hdr.get('INSTRUME', 'WSRT / Apertif'),
        survey             = hdr.get('OBJECT', dataset_name),
        dataset            = dataset_name,
        provider           = 'ASTRON',
        observation_id     = dataset_name,
        observatory_name   = 'WSRT',
        latitude_deg       = WSRT_LAT_DEG,
        longitude_deg      = WSRT_LON_DEG,
        elevation_m        = WSRT_ELEV_M,
        frequency_axis_hz  = frequency_axis_hz,
        center_frequency_hz= center_hz,
        receiver           = hdr.get('INSTRUME', 'WSRT / Apertif'),
        pointing_source    = 'measured',
        notes              = f"WSRT/Apertif FITS cube {fits_path.name}. "
                             f"Integrated spectrum (nanmean over spatial axes). "
                             f"Parsed by wsrt_to_fro.py.",
        overwrite          = False
    )

    fro.append_raw_spectrum(
        filename           = str(fro_path),
        raw_spectrum       = raw_spectrum,
        utc_time           = utc_dt,
        integration_time_s = float(hdr.get('EXPTIME', float('nan'))),
        az_deg             = float('nan'),
        elev_deg           = float('nan'),
        ra_deg             = ra_deg,
        dec_deg            = dec_deg,
        glon_deg           = glon_deg,
        glat_deg           = glat_deg,
    )

    print(f"FRO scritto: {fro_path}")
    return fro_path


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  WSRT / Apertif → HALO v17 parser")
    print("=" * 60)

    if len(sys.argv) > 1:
        input_path = Path(sys.argv[1])
    else:
        input_path = select_input()

    if not input_path.exists():
        print(f"Percorso non trovato: {input_path}")
        sys.exit(1)

    mode = detect_mode(input_path)

    if mode == 'ms':
        parse_ms(input_path)
    elif mode == 'uvfits':
        parse_uvfits(input_path)
    else:
        parse_fits(input_path)

    print("\nDone.")


if __name__ == '__main__':
    main()
