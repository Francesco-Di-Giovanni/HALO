#!/usr/bin/env python3
"""
meerkat_to_fro.py — MeerKAT to FRO/HALO v17 parser
Supports two input modes:
  - MeasurementSet (.ms directory): extracts autocorrelations via casatools
  - FITS cube (.fits): extracts spectra by position via astropy

Autocorrelations (ANTENNA1 == ANTENNA2) in a MeasurementSet are equivalent
to single-dish total-power spectra for each antenna and integration.

Usage:
    python3 meerkat_to_fro.py

Requires:
  - astropy, numpy, h5py (standard FRO environment)
  - casatools (conda env 'casa_env') for .ms input mode
  - fro_format_v17.py in ~/FRO/FRO_System/

Output:
  ~/FRO/Parsers/MeerKAT/Surveys/<dataset>_<YYYYMMDD_HHMMSS>/FRO/<name>.fro

Author: Francesco Di Giovanni / FRO-HALO project
"""

import os
import sys
import subprocess
import numpy as np
from datetime import datetime, timezone
from pathlib import Path

# ── FRO system path ──────────────────────────────────────────────────────────
FRO_SYSTEM = Path.home() / 'FRO' / 'FRO_System'
sys.path.insert(0, str(FRO_SYSTEM))
import fro_format_v17 as fro

# ── Output base directory ────────────────────────────────────────────────────
SURVEYS_DIR = Path.home() / 'FRO' / 'Parsers' / 'MeerKAT' / 'Surveys'

# ── MeerKAT observatory coordinates (core array, Karoo, South Africa) ────────
MEERKAT_LAT_DEG  =  -30.7130        # deg N
MEERKAT_LON_DEG  =   21.4430        # deg E
MEERKAT_ELEV_M   =  1035.0          # m a.s.l.

# ── Sentinel for unknown epoch (mosaic/survey products) ──────────────────────
UNKNOWN_EPOCH = fro.UNKNOWN_EPOCH


# ════════════════════════════════════════════════════════════════════════════
# Input selection
# ════════════════════════════════════════════════════════════════════════════

def select_input():
    """Prompt user to select input file or directory via zenity."""
    try:
        result = subprocess.run(
            ['zenity', '--file-selection',
             '--title=Seleziona MeasurementSet (.ms) o FITS cube',
             '--filename=' + str(Path.home() / 'Downloads') + '/',
             '--file-filter=MeasurementSet (directory) | *.ms',
             '--file-filter=FITS cube | *.fits *.FITS'],
            capture_output=True, text=True
        )
        path = result.stdout.strip()
        if not path:
            print("Nessun file selezionato. Uscita.")
            sys.exit(0)
        return Path(path)
    except FileNotFoundError:
        # Fallback: command line
        path = input("Percorso MeasurementSet o FITS cube: ").strip()
        if not path:
            sys.exit(0)
        return Path(path)


def detect_mode(input_path):
    """Return 'ms' or 'fits' based on input path."""
    p = str(input_path).lower()
    if input_path.is_dir() and (p.endswith('.ms') or (input_path / 'table.info').exists()):
        return 'ms'
    if input_path.is_file() and (p.endswith('.fits') or p.endswith('.fit')):
        return 'fits'
    # Directory without .ms extension: assume MeasurementSet
    if input_path.is_dir():
        return 'ms'
    print(f"Tipo input non riconosciuto: {input_path}")
    sys.exit(1)


# ════════════════════════════════════════════════════════════════════════════
# Mode A: MeasurementSet via casatools
# ════════════════════════════════════════════════════════════════════════════

def parse_ms(ms_path):
    """
    Extract autocorrelation spectra from a MeerKAT MeasurementSet.
    Autocorrelations (ANTENNA1 == ANTENNA2) are single-dish total-power spectra.
    Uses casatools from the casa_env conda environment.
    """
    print(f"Modalità: MeasurementSet")
    print(f"Input: {ms_path}")

    # casatools is available in casa_env; import at runtime
    try:
        from casatools import ms as mstool, table as tbtool
    except ImportError:
        print("\nErrore: casatools non trovato.")
        print("Attivare l'ambiente conda: conda activate casa_env")
        print("Poi rieseguire: python3 meerkat_to_fro.py")
        sys.exit(1)

    ms = mstool()
    tb = tbtool()

    ms.open(str(ms_path))
    ms_info = ms.getscansummary()
    spw_info = ms.getspectralwindowinfo()
    ms.close()

    # Read antenna positions and names
    tb.open(str(ms_path) + '/ANTENNA')
    antenna_names = list(tb.getcol('NAME'))
    antenna_pos   = tb.getcol('POSITION')   # ITRF XYZ, shape (3, N_ant)
    tb.close()

    # Read spectral window(s)
    tb.open(str(ms_path) + '/SPECTRAL_WINDOW')
    chan_freq  = tb.getcol('CHAN_FREQ')      # shape (N_chan, N_spw)
    chan_width = tb.getcol('CHAN_WIDTH')     # shape (N_chan, N_spw)
    ref_freq   = tb.getcol('REF_FREQUENCY') # shape (N_spw,)
    tb.close()

    # Read field/source info
    tb.open(str(ms_path) + '/FIELD')
    field_names  = list(tb.getcol('NAME'))
    phase_dir    = tb.getcol('PHASE_DIR')   # (2, 1, N_field) rad (RA, Dec J2000)
    tb.close()

    # Read main table — autocorrelations only
    tb.open(str(ms_path))
    ant1     = tb.getcol('ANTENNA1')
    ant2     = tb.getcol('ANTENNA2')
    data     = tb.getcol('DATA')            # (N_corr, N_chan, N_row) complex
    time_col = tb.getcol('TIME')            # MJD seconds
    field_id = tb.getcol('FIELD_ID')
    spw_id   = tb.getcol('DATA_DESC_ID')   # proxy for spw
    tb.close()

    # Filter autocorrelations
    auto_mask = (ant1 == ant2)
    if not np.any(auto_mask):
        print("Nessuna autocorrelazione trovata nel MeasurementSet.")
        sys.exit(1)

    data_auto     = data[:, :, auto_mask]       # (N_corr, N_chan, N_auto)
    time_auto     = time_col[auto_mask]          # MJD seconds
    ant1_auto     = ant1[auto_mask]
    field_auto    = field_id[auto_mask]
    spw_auto      = spw_id[auto_mask]

    n_auto = data_auto.shape[2]
    print(f"Autocorrelazioni trovate: {n_auto}")

    # Determine dataset name from MS path
    dataset_name = ms_path.stem  # e.g. "1234567890.ms" → "1234567890"

    # Epoch directory
    epoch_str = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    epoch_dir = SURVEYS_DIR / f"{dataset_name}_{epoch_str}"
    fro_dir   = epoch_dir / 'FRO'
    raw_dir   = epoch_dir / 'Raw'
    fro_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    # Write source reference
    with open(raw_dir / 'source.txt', 'w') as f:
        f.write(f"MeasurementSet: {ms_path}\n")
        f.write(f"Antennas: {', '.join(antenna_names)}\n")

    # Use first spw for frequency axis (extend to multi-spw if needed)
    spw_idx   = int(spw_auto[0])
    frequency_axis_hz   = chan_freq[:, spw_idx].astype(np.float64)   # Hz
    n_chan    = len(frequency_axis_hz)
    center_hz = float(ref_freq[spw_idx])

    # Integration time: difference between consecutive timestamps (same antenna)
    unique_times = np.unique(time_auto)
    if len(unique_times) > 1:
        integration_s = float(np.median(np.diff(unique_times)))
    else:
        integration_s = float('nan')

    # Field 0 coordinates (primary field)
    ra_rad  = float(phase_dir[0, 0, 0])
    dec_rad = float(phase_dir[1, 0, 0])
    ra_deg  = np.degrees(ra_rad)
    dec_deg = np.degrees(dec_rad)

    # Convert to Galactic
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    gal   = coord.galactic
    glon_deg = float(gal.l.deg)
    glat_deg = float(gal.b.deg)

    # Output .fro path
    fro_path = fro_dir / f"{dataset_name}_{epoch_str}.fro"

    # Create FRO structure
    fro.create_fro_v17_structure(
        filename           = str(fro_path),
        facility        = 'MeerKAT',
        instrument      = 'MeerKAT L-band',
        survey          = dataset_name,
        dataset         = dataset_name,
        provider        = 'SARAO',
        observation_id  = dataset_name,
        observatory_name= 'MeerKAT',
        lat_deg         = MEERKAT_LAT_DEG,
        lon_deg         = MEERKAT_LON_DEG,
        elev_m          = MEERKAT_ELEV_M,
        frequency_axis_hz         = frequency_axis_hz,
        center_frequency_axis_hz  = center_hz,
        receiver        = 'L-band (856–1712 MHz)',
        integration_time_s = n_auto,
        pointing_source = 'measured',
        notes           = f"MeerKAT autocorrelations from {ms_path.name}. "
                          f"Parsed by meerkat_to_fro.py. "
                          f"ANTENNA1==ANTENNA2 rows only.",
        overwrite       = False
    )

    # Append spectra
    # Use Stokes I approximation: mean of XX and YY (corr 0 and 3 if present)
    n_corr = data_auto.shape[0]
    if n_corr >= 4:
        stokes_i = 0.5 * (np.abs(data_auto[0, :, :]) + np.abs(data_auto[3, :, :]))
    else:
        stokes_i = np.abs(data_auto[0, :, :])   # fallback: single corr

    for i in range(n_auto):
        # MJD seconds → UTC datetime
        mjd_s   = time_auto[i]
        mjd_day = mjd_s / 86400.0
        from astropy.time import Time
        t = Time(mjd_day, format='mjd', scale='utc')
        utc_dt  = t.to_datetime(timezone.utc)

        raw_spectrum = stokes_i[:, i].astype(np.float32)

        fro.append_raw_spectrum(
            filename           = str(fro_path),
            spectrum         = spectrum,
            utc_time         = utc_dt,
            integration_time_s = integration_s,
            azimuth_deg      = float('nan'),
            elevation_deg    = float('nan'),
            ra_deg           = ra_deg,
            dec_deg          = dec_deg,
            glon_deg         = glon_deg,
            glat_deg         = glat_deg,
        )

        if (i + 1) % 100 == 0 or i == n_auto - 1:
            print(f"  Spettro {i+1}/{n_auto}", end='\r')

    print(f"\nFRO scritto: {fro_path}")
    return fro_path


# ════════════════════════════════════════════════════════════════════════════
# Mode B: FITS cube
# ════════════════════════════════════════════════════════════════════════════

def parse_fits(fits_path):
    """
    Extract spectra from a MeerKAT FITS cube (reduced survey product).
    Cube axes expected: (STOKES, FREQ, DEC, RA) or (FREQ, DEC, RA).
    Extracts integrated spectrum over all spatial pixels.
    """
    print(f"Modalità: FITS cube")
    print(f"Input: {fits_path}")

    from astropy.io import fits as astrofits
    from astropy.coordinates import SkyCoord, Galactic
    from astropy.wcs import WCS
    import astropy.units as u
    import warnings
    from astropy.utils.exceptions import AstropyWarning

    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        with astrofits.open(str(fits_path)) as hdul:
            hdr  = hdul[0].header
            data = hdul[0].data  # float array

    # WCS for frequency and spatial axes
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', AstropyWarning)
        wcs = WCS(hdr)

    # Determine frequency axis index (CTYPE = 'FREQ')
    n_axes = data.ndim
    freq_axis = None
    for i in range(1, n_axes + 1):
        ct = hdr.get(f'CTYPE{i}', '')
        if 'FREQ' in ct.upper():
            freq_axis = i - 1   # 0-based, FITS order
            break

    if freq_axis is None:
        print("Asse FREQ non trovato nell'header FITS.")
        sys.exit(1)

    n_freq = data.shape[n_axes - 1 - freq_axis]

    # Build frequency array
    crpix = hdr.get(f'CRPIX{freq_axis+1}', 1.0)
    crval = hdr.get(f'CRVAL{freq_axis+1}', 0.0)
    cdelt = hdr.get(f'CDELT{freq_axis+1}', 1.0)
    frequency_axis_hz = (crval + (np.arange(n_freq) - (crpix - 1)) * cdelt).astype(np.float64)
    center_hz = float(crval)

    # Integration / epoch
    date_obs = hdr.get('DATE-OBS', None)
    if date_obs:
        try:
            from astropy.time import Time
            t = Time(date_obs, format='isot', scale='utc')
            utc_dt = t.to_datetime(timezone.utc)
        except Exception:
            utc_dt = UNKNOWN_EPOCH
    else:
        utc_dt = UNKNOWN_EPOCH

    # Pointing: phase centre from WCS
    try:
        ra_deg  = float(hdr.get('CRVAL1', float('nan')))
        dec_deg = float(hdr.get('CRVAL2', float('nan')))
        coord   = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
        gal     = coord.galactic
        glon_deg = float(gal.l.deg)
        glat_deg = float(gal.b.deg)
    except Exception:
        ra_deg = dec_deg = glon_deg = glat_deg = float('nan')

    # Integrated spectrum: mean over all spatial axes, ignoring NaN
    # Collapse all axes except frequency
    # Move freq axis to last position, then mean over all others
    freq_numpy_axis = n_axes - 1 - freq_axis  # numpy axis index
    axes_to_mean = tuple(i for i in range(n_axes) if i != freq_numpy_axis)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        raw_spectrum = np.nanmean(data, axis=axes_to_mean).astype(np.float32)

    # Receiver from header
    receiver = hdr.get('INSTRUME', 'MeerKAT L-band')

    # Dataset name
    dataset_name = fits_path.stem

    # Epoch directory
    epoch_str = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    epoch_dir = SURVEYS_DIR / f"{dataset_name}_{epoch_str}"
    fro_dir   = epoch_dir / 'FRO'
    raw_dir   = epoch_dir / 'Raw'
    fro_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    with open(raw_dir / 'source.txt', 'w') as f:
        f.write(f"FITS cube: {fits_path}\n")

    fro_path = fro_dir / f"{dataset_name}_{epoch_str}.fro"

    fro.create_fro_v17_structure(
        filename           = str(fro_path),
        facility        = 'MeerKAT',
        instrument      = receiver,
        survey          = hdr.get('SURVEY', dataset_name),
        dataset         = dataset_name,
        provider        = hdr.get('ORIGIN', 'SARAO'),
        observation_id  = dataset_name,
        observatory_name= 'MeerKAT',
        lat_deg         = MEERKAT_LAT_DEG,
        lon_deg         = MEERKAT_LON_DEG,
        elev_m          = MEERKAT_ELEV_M,
        frequency_axis_hz         = frequency_axis_hz,
        center_frequency_axis_hz  = center_hz,
        receiver        = receiver,
        integration_time_s = 1,
        pointing_source = 'measured',
        notes           = f"MeerKAT FITS cube {fits_path.name}. "
                          f"Integrated spectrum (nanmean over spatial axes). "
                          f"Parsed by meerkat_to_fro.py.",
        overwrite       = False
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
    print("  MeerKAT → FRO/HALO v17 parser")
    print("=" * 60)

    input_path = select_input()

    if not input_path.exists():
        print(f"Percorso non trovato: {input_path}")
        sys.exit(1)

    mode = detect_mode(input_path)

    if mode == 'ms':
        fro_path = parse_ms(input_path)
    else:
        fro_path = parse_fits(input_path)

    print("\nDone.")
    print(f"Output: {fro_path}")


if __name__ == '__main__':
    main()
