#!/usr/bin/env python3
"""
things_to_fro.py — THINGS (The HI Nearby Galaxy Survey) to HALO v17 parser

THINGS is a VLA HI survey of 34 nearby galaxies (Walter et al. 2008).
Data products are FITS cubes with axes (RA, Dec, Velocity) or (RA, Dec, Freq).
Natural-weighted and robust-weighted cubes are available.

This parser extracts an integrated HI spectrum from a THINGS FITS cube
by averaging over all spatial pixels (nanmean), preserving the spectral axis.

Usage:
    python3 things_to_fro.py [input.fits]

Requires:
  - astropy, numpy, h5py (standard HALO environment)
  - fro_format_v17.py in ~/FRO/FRO_System/

Output:
  ~/FRO/Parsers/THINGS/Surveys/<dataset>_<YYYYMMDD_HHMMSS>/FRO/<name>.fro

Data source: https://www.mpia.de/THINGS/Data.html
Reference: Walter et al. 2008, AJ, 136, 2563

Project: HALO (Hydrogen Atomic Line Observatory)
Author: Francesco Di Giovanni (Bolzano, Italy)
Created: 2026-08-21 (Parte 67)
"""

import os
import sys
import numpy as np
from datetime import datetime, timezone
from pathlib import Path

# ── HALO system path ─────────────────────────────────────────────────────────
FRO_SYSTEM = Path.home() / 'FRO' / 'FRO_System'
sys.path.insert(0, str(FRO_SYSTEM))
import fro_format_v17 as fro

# ── Output base directory ────────────────────────────────────────────────────
SURVEYS_DIR = Path.home() / 'FRO' / 'Parsers' / 'THINGS' / 'Surveys'

# ── VLA coordinates (THINGS observed with VLA, B/C/D configurations) ────────
# VLA is located at Socorro, New Mexico, USA
VLA_LAT_DEG =  34.0784
VLA_LON_DEG = -107.6184
VLA_ELEV_M  = 2124.0

# ── HI rest frequency ────────────────────────────────────────────────────────
HI_REST_HZ = 1420405751.768

# ── Sentinel for unknown epoch ───────────────────────────────────────────────
UNKNOWN_EPOCH = fro.UNKNOWN_EPOCH


def build_freq_axis(hdr, n_spec, spec_axis_fits):
    """
    Build frequency axis in Hz from FITS header.
    Handles both FREQ and VRAD/VOPT axes (converting velocity to frequency).
    spec_axis_fits: 1-based FITS axis index of the spectral axis.
    """
    ctype = hdr.get(f'CTYPE{spec_axis_fits}', '').upper().strip()
    crpix = float(hdr.get(f'CRPIX{spec_axis_fits}', 1.0))
    crval = float(hdr.get(f'CRVAL{spec_axis_fits}', 0.0))
    cdelt = float(hdr.get(f'CDELT{spec_axis_fits}', 1.0))

    if 'FREQ' in ctype:
        # Direct frequency axis (Hz)
        frequency_axis_hz = crval + (np.arange(n_spec) - (crpix - 1)) * cdelt
        center_hz = float(crval)
    elif 'VRAD' in ctype or 'VOPT' in ctype or 'VELO' in ctype or 'FELO' in ctype:
        # Velocity axis (m/s) — convert to frequency using HI rest frequency
        # Radio convention: f = f0 * (1 - v/c)
        c_mps = 299792458.0
        vel_mps = crval + (np.arange(n_spec) - (crpix - 1)) * cdelt
        frequency_axis_hz = HI_REST_HZ * (1.0 - vel_mps / c_mps)
        center_hz = float(HI_REST_HZ * (1.0 - crval / c_mps))
        print(f"Asse spettrale: {ctype} — convertito in frequenza via HI rest {HI_REST_HZ/1e6:.6f} MHz")
    else:
        print(f"WARNING: CTYPE spettrale non riconosciuto: '{ctype}'. Asse in unità originali.")
        frequency_axis_hz = crval + (np.arange(n_spec) - (crpix - 1)) * cdelt
        center_hz = float(crval)

    return frequency_axis_hz.astype(np.float64), center_hz


def parse_fits(fits_path):
    """Extract integrated HI spectrum from a THINGS FITS cube."""
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
            data = hdul[0].data  # shape e.g. (1, n_vel, n_dec, n_ra) or (n_vel, n_dec, n_ra)

    # Squeeze degenerate Stokes axis if present
    if data.ndim == 4:
        data = data.squeeze(axis=0)  # remove Stokes axis (always I for THINGS)

    n_axes = data.ndim  # should be 3: (vel/freq, dec, ra)

    # Find spectral axis (FREQ, VRAD, VOPT, VELO, FELO)
    spec_axis_fits = None
    for i in range(1, hdr.get('NAXIS', 0) + 1):
        ct = hdr.get(f'CTYPE{i}', '').upper()
        if any(k in ct for k in ('FREQ', 'VRAD', 'VOPT', 'VELO', 'FELO')):
            spec_axis_fits = i
            break

    if spec_axis_fits is None:
        print("Asse spettrale non trovato nell'header FITS.")
        sys.exit(1)

    # In numpy, FITS axis order is reversed: axis 0 = NAXIS_last
    # For a 3D cube (NAXIS1=RA, NAXIS2=Dec, NAXIS3=Vel):
    # numpy shape = (n_vel, n_dec, n_ra), spectral numpy axis = 0
    spec_numpy_axis = hdr.get('NAXIS', 0) - spec_axis_fits
    # Adjust for squeezed Stokes axis
    if data.ndim == 3 and hdr.get('NAXIS', 0) == 4:
        spec_numpy_axis = max(0, spec_numpy_axis - 1)

    n_spec = data.shape[spec_numpy_axis]
    frequency_axis_hz, center_hz = build_freq_axis(hdr, n_spec, spec_axis_fits)

    # Integrated spectrum: nanmean over all spatial axes
    spatial_axes = tuple(i for i in range(data.ndim) if i != spec_numpy_axis)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        raw_spectrum = np.nanmean(data, axis=spatial_axes).astype(np.float32)

    # Ensure spectrum aligns with ascending frequency axis
    if len(frequency_axis_hz) > 1 and frequency_axis_hz[0] > frequency_axis_hz[-1]:
        frequency_axis_hz  = frequency_axis_hz[::-1]
        spectrum = spectrum[::-1]

    # Pointing: phase centre
    try:
        ra_deg  = float(hdr.get('CRVAL1', float('nan')))
        dec_deg = float(hdr.get('CRVAL2', float('nan')))
        coord   = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
        gal     = coord.galactic
        glon_deg = float(gal.l.deg)
        glat_deg = float(gal.b.deg)
    except Exception:
        ra_deg = dec_deg = glon_deg = glat_deg = float('nan')

    # Epoch
    date_obs = hdr.get('DATE-OBS', None)
    if date_obs:
        try:
            from astropy.time import Time
            utc_dt = Time(date_obs, format='isot', scale='utc').to_datetime(timezone.utc)
        except Exception:
            utc_dt = UNKNOWN_EPOCH
    else:
        utc_dt = UNKNOWN_EPOCH

    # Source name from header
    object_name = hdr.get('OBJECT', fits_path.stem).strip()
    dataset_name = fits_path.stem
    epoch_str    = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S')
    epoch_dir    = SURVEYS_DIR / f"{dataset_name}_{epoch_str}"
    fro_dir      = epoch_dir / 'FRO'
    raw_dir      = epoch_dir / 'Raw'
    fro_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)

    with open(raw_dir / 'source.txt', 'w') as f:
        f.write(f"THINGS FITS cube: {fits_path}\n")
        f.write(f"Object: {object_name}\n")
        f.write(f"Reference: Walter et al. 2008, AJ, 136, 2563\n")

    fro_path = fro_dir / f"{dataset_name}_{epoch_str}.fro"

    bunit = hdr.get('BUNIT', 'unknown').strip()

    fro.create_fro_v17_structure(
        filename           = str(fro_path),
        facility           = 'VLA',
        instrument         = 'VLA L-band',
        survey             = 'THINGS',
        dataset            = dataset_name,
        provider           = 'NRAO / MPIA',
        observation_id     = dataset_name,
        observatory_name   = 'VLA',
        latitude_deg       = VLA_LAT_DEG,
        longitude_deg      = VLA_LON_DEG,
        elevation_m        = VLA_ELEV_M,
        frequency_axis_hz  = frequency_axis_hz,
        center_frequency_hz= center_hz,
        receiver           = 'VLA L-band (HI 21cm)',
        pointing_source    = 'measured',
        notes              = f"THINGS HI cube {fits_path.name}. "
                             f"Object: {object_name}. "
                             f"Integrated spectrum (nanmean over spatial axes). "
                             f"Units: {bunit}. "
                             f"Walter et al. 2008, AJ, 136, 2563. "
                             f"Parsed by things_to_fro.py.",
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

    print("=" * 60)
    print("  THINGS → HALO v17 conversion complete")
    print("=" * 60)
    print(f"Input       : {fits_path}")
    print(f"Object      : {object_name}")
    print(f"Channels    : {n_spec}")
    print(f"Freq range  : {frequency_axis_hz[0]/1e6:.3f} – {frequency_axis_hz[-1]/1e6:.3f} MHz")
    print(f"Units       : {bunit}")
    print(f"FRO scritto : {fro_path}")
    print("=" * 60)
    return fro_path


def main():
    print("=" * 60)
    print("  THINGS → HALO v17 parser")
    print("=" * 60)

    if len(sys.argv) > 1:
        fits_path = Path(sys.argv[1])
    else:
        import subprocess
        try:
            result = subprocess.run(
                ['zenity', '--file-selection',
                 '--title=Seleziona THINGS FITS cube',
                 '--filename=' + str(Path.home() / 'Downloads') + '/',
                 '--file-filter=FITS files | *.fits *.fit'],
                capture_output=True, text=True
            )
            path = result.stdout.strip()
            if not path:
                print("Nessun file selezionato. Uscita.")
                sys.exit(0)
            fits_path = Path(path)
        except FileNotFoundError:
            path = input("Percorso FITS cube: ").strip()
            if not path:
                sys.exit(0)
            fits_path = Path(path)

    if not fits_path.exists():
        print(f"File non trovato: {fits_path}")
        sys.exit(1)

    parse_fits(fits_path)
    print("\nDone.")


if __name__ == '__main__':
    main()
