#!/usr/bin/env python3
# =============================================================================
# fast_fits_to_fro.py
#
# Parser: FAST raw FITS file -> FRO format v17 (HDF5)
#
# Project: FRO (Francesco Radio Observatory) / HALO
# Author: Francesco Di Giovanni (Bolzano, Italy)
# Created: 2026-07-13 (Parte 35)
#
# Reads a FAST L-band 19-beam raw FITS file (one file per beam, chunked)
# and converts it to FRO v17 HDF5 using fro_format_v17.py.
#
# FAST raw FITS structure (verified on M33_OTF_1_MultiBeamOTF-M01_W_0001.fits):
#   - PRIMARY: empty (4 cards only)
#   - SINGLE DISH BinTable: 615R x 21C
#     KEY COLUMNS:
#       DATE-OBS  : ISO timestamp string with 'Z' suffix (UTC)
#       UTOBS     : MJD (double)
#       OBJ_RA    : RA of object (deg) -- ZERO in raw files (no pointing)
#       OBJ_DEC   : Dec of object (deg) -- ZERO in raw files (no pointing)
#       OBSTYPE   : 'ON'/'OFF'/'CAL' -> flag_cal
#       QUALITY   : boolean -> flag_track
#       TSYS      : system temperature (K) -- zero in raw files
#       EXPOSURE  : integration time per dump (s) ~0.5s
#       NCHAN     : number of channels (65536)
#       FREQ      : frequency of first channel (MHz)
#       CHAN_BW   : channel bandwidth (MHz)
#       PRESSURE  : atmospheric pressure -- zero in raw files
#       TAMBIENT  : ambient temperature -- zero in raw files
#       WINDSPD   : wind speed
#       WINDDIR   : wind direction
#       DATA      : (65536, 4) float32 array
#                   columns: XX, YY, Re(XY), Im(XY)
#                   (CTYPE2=STOKES, CRVAL2=-5.0 -> XX=-5, YY=-6, XY=-7, YX=-8)
#                   TDIM21=(4,65536) Fortran order -> numpy shape (65536,4)
#
# POINTING NOTE (v1):
#   OBJ_RA/OBJ_DEC are zero in all rows of the raw FITS. Celestial
#   coordinates require the KY .xlsx file (feed cabin mechanical positions)
#   and the FAST geometric model (hifast.radec). All pointing fields are
#   stored as NaN in this v1 import. Future v2 will implement coordinate
#   reconstruction when the geometric model is available.
#
# OUTPUT: one .fro file per input FITS file (one beam, one chunk).
#         XX and YY stored as separate spectra (two rows per dump).
#         Re(XY) and Im(XY) stored as cross spectra (include_cross=False
#         in v1 -- cross pol support can be added in v2).
# =============================================================================

import os
import sys
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
from scipy.interpolate import interp1d
import astropy.units as u

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fro_format_v17 import (
    create_fro_v17_structure,
    append_raw_spectrum,
)

# ---------------------------------------------------------------------------
# FAST site constants (Pingtang, Guizhou)
# ---------------------------------------------------------------------------
FAST_SITE_LAT_DEG = 25.6525
FAST_SITE_LON_DEG = 106.8569
FAST_SITE_ELEV_M  = 1110.0
FAST_TELESCOPE    = "FAST"
FAST_FRONTEND     = "L-band 19-beam"
# ---------------------------------------------------------------------------
# KY coordinate table loader
# ---------------------------------------------------------------------------
FAST_LOCATION = EarthLocation.from_geodetic(
    FAST_SITE_LON_DEG * u.deg,
    FAST_SITE_LAT_DEG * u.deg,
    height=FAST_SITE_ELEV_M * u.m,
)

def load_ky_radec(ky_table_path, beam=1):
    """Load KY RA/Dec table and return interpolators for given beam.
    Table format: MJD RA1 Dec1 RA2 Dec2 ... RA19 Dec19
    """
    d = np.loadtxt(ky_table_path)
    mjd = d[:, 0]
    ra  = d[:, (beam - 1) * 2 + 1]
    dec = d[:, (beam - 1) * 2 + 2]
    ra_interp  = interp1d(mjd, ra,  kind='linear',
                          bounds_error=False, fill_value=float('nan'))
    dec_interp = interp1d(mjd, dec, kind='linear',
                          bounds_error=False, fill_value=float('nan'))
    return ra_interp, dec_interp

# FITS STOKES polarization codes (CRVAL2 axis)
POL_CODES = {
    -1: "RR", -2: "LL", -3: "RL", -4: "LR",
    -5: "XX", -6: "YY", -7: "XY", -8: "YX",
}


def extract_beam_number(fits_path):
    """Extract beam number from FAST file naming convention.
    Expected: ..._MultiBeamOTF-M{beam:02d}_W_{chunk:04d}.fits
    """
    basename = os.path.basename(fits_path)
    try:
        part = basename.split("MultiBeamOTF-M")[1]
        return int(part.split("_")[0])
    except (IndexError, ValueError):
        return -1


def extract_chunk_number(fits_path):
    """Extract chunk number from FAST file naming convention."""
    basename = os.path.basename(fits_path)
    try:
        part = basename.split("_W_")[1]
        return int(part.split(".")[0])
    except (IndexError, ValueError):
        return -1


def parse_date_obs(date_obs_str):
    """Convert DATE-OBS string to UTC datetime.
    FAST format: '2021-07-30T21:34:07.959Z' (ISO with Z suffix).
    """
    s = str(date_obs_str).strip().rstrip("Z")
    dt = datetime.fromisoformat(s)
    return dt.replace(tzinfo=timezone.utc)


def build_freq_axis(freq_mhz, chan_bw_mhz, nchan):
    """Build sky frequency axis in Hz.
    FREQ = frequency of first channel (MHz), CHAN_BW = channel width (MHz).
    Axis is ascending (CHAN_BW > 0 verified on reference file).
    """
    return (freq_mhz + np.arange(nchan) * chan_bw_mhz) * 1e6


def convert_fast_fits(fits_path, out_dir=None, ky_table_path=None):
    """Convert one FAST raw FITS file (one beam, one chunk) to .fro v17."""
    fits_path = os.path.abspath(fits_path)
    if out_dir is None:
        out_dir = os.path.dirname(fits_path)

    beam_number  = extract_beam_number(fits_path)
    chunk_number = extract_chunk_number(fits_path)
    # --- KY coordinate interpolators ---
    if ky_table_path is not None and os.path.exists(ky_table_path):
        ra_interp, dec_interp = load_ky_radec(ky_table_path, beam=beam_number)
        has_ky = True
    else:
        ra_interp = dec_interp = None
        has_ky = False

    hdul = fits.open(fits_path)

    # Find SINGLE DISH extension
    ext = None
    for i, hdu in enumerate(hdul):
        if hasattr(hdu, "name") and "SINGLE DISH" in str(hdu.name).upper():
            ext = hdu
            break
    if ext is None:
        hdul.close()
        raise ValueError("No 'SINGLE DISH' extension found in FITS file.")

    data = ext.data
    hdr  = ext.header
    n_rows = len(data)

    # --- Frequency axis (from first row, verified stable) ---
    freq_mhz   = float(data[0]["FREQ"])
    chan_bw_mhz = float(data[0]["CHAN_BW"])
    nchan      = int(data[0]["NCHAN"])
    freq_axis_hz = build_freq_axis(freq_mhz, chan_bw_mhz, nchan)

    # --- Polarization labels from STOKES WCS ---
    crval2 = float(hdr.get("CRVAL2", -5.0))  # -5 = XX
    pol0 = POL_CODES.get(int(crval2),     f"POL{int(crval2)}")
    pol1 = POL_CODES.get(int(crval2) - 1, f"POL{int(crval2)-1}")

    # --- Global metadata ---
    target   = hdr.get("OBJECT",   f"beam{beam_number:02d}")
    projid   = hdr.get("PROJID",   "unknown")
    observer = hdr.get("OBSERVER", "unknown")
    exposure_s = float(data[0]["EXPOSURE"])
    bandwidth_hz = nchan * chan_bw_mhz * 1e6

    # --- Output file ---
    base     = os.path.splitext(os.path.basename(fits_path))[0]
    out_path = os.path.join(out_dir, base + "_FRO_v17.h5")
    observation_id = (f"{FAST_TELESCOPE}_{projid}"
                      f"_beam{beam_number:02d}_chunk{chunk_number:04d}")

    notes = (
        f"Converted from FAST raw FITS file '{os.path.basename(fits_path)}'. "
        f"Beam {beam_number}, chunk {chunk_number}, {n_rows} dumps. "
        f"Observer: {observer}. Project: {projid}. "
        "POINTING (v1): OBJ_RA/OBJ_DEC are zero in all rows of the raw "
        "FITS file. Celestial coordinates require the KY .xlsx file "
        "(feed cabin mechanical positions, 3D TRP/SDP coordinate systems) "
        "and the FAST geometric model (hifast.radec) for conversion to "
        "RA/DEC. All pointing fields stored as NaN in this v1 import. "
        "DATA array: shape (65536, 4), columns XX/YY/Re(XY)/Im(XY) "
        "(CTYPE2=STOKES, CRVAL2=-5, TDIM=(4,65536) Fortran order). "
        "Frequency axis: FREQ + i*CHAN_BW (MHz), verified linear ascending. "
        "OBSTYPE column: ON/OFF/CAL -> stored in flag_cal "
        "(0=ON, 1=CAL, -1=OFF/unknown). "
        "PRESSURE/TAMBIENT/WINDSPD/WINDDIR: zero in this raw file "
        "(not populated at raw level), stored as NaN."
    )

    temperature_calibration_note = (
        "Raw uncalibrated counts (float32). TSYS=0 in raw FAST files. "
        "Calibration requires noise diode injection phases (CAL rows in "
        "OBSTYPE) and frequency-dependent Tcal from FAST calibration tables."
    )

    create_fro_v17_structure(
        out_path,
        frequency_axis_hz=freq_axis_hz,
        facility=FAST_TELESCOPE,
        instrument=FAST_FRONTEND,
        observation_id=observation_id,
        survey="",
        dataset="",
        provider="FAST public data (HiFAST sample dataset)",
        import_date=datetime.now(timezone.utc).isoformat(),
        observatory_name=FAST_TELESCOPE,
        latitude_deg=FAST_SITE_LAT_DEG,
        longitude_deg=FAST_SITE_LON_DEG,
        elevation_m=FAST_SITE_ELEV_M,
        receiver=FAST_FRONTEND,
        antenna=FAST_TELESCOPE,
        lna="n/a",
        center_frequency_hz=(freq_mhz + nchan * chan_bw_mhz / 2.0) * 1e6,
        sample_rate_hz=float("nan"),
        bandwidth_hz=bandwidth_hz,
        integration_time_s=exposure_s,
        notes=notes,
        temperature_calibration_note=temperature_calibration_note,
        processed_by="fast_fits_to_fro.py",
        pointing_source="NaN - see Notes",
        target=target,
        vlsr_mps=float("nan"),
        scan_id=chunk_number,
        subscan_id=beam_number,
        subscan_type="OTF",
        signal="unknown",
        schedule_name=projid,
        cal_polarizations=[pol0, pol1],
        cal_mark_temp_k=[float("nan"), float("nan")],
        attenuation_db=[float("nan"), float("nan")],
        include_cross=False,
    )

    # --- Per-dump loop ---
    print(f"Converting {n_rows} dumps, beam {beam_number}, chunk {chunk_number}...")
    for i in range(n_rows):
        row  = data[i]
        arr  = np.asarray(row["DATA"], dtype=np.float32)  # shape (65536, 4)
        xx   = arr[:, 0]
        yy   = arr[:, 1]

        utc_dt = parse_date_obs(row["DATE-OBS"])

        # flag_cal: 0=ON (science), 1=CAL (calibration), -1=OFF/unknown
        obstype = str(row["OBSTYPE"]).strip().upper()
        if obstype == "ON":
            flag_cal = 0
        elif obstype == "CAL":
            flag_cal = 1
        else:
            flag_cal = -1

        # flag_track: from QUALITY boolean
        flag_track = 1 if bool(row["QUALITY"]) else 0

        # Environment: all zero in raw file -> NaN
        tambient = float(row["TAMBIENT"])
        pressure = float(row["PRESSURE"])
        temperature_c   = tambient - 273.15 if tambient > 0 else float("nan")
        pressure_mbar   = pressure if pressure > 0 else float("nan")

        # Two rows per dump (one per polarization)
        for spec, pol in ((xx, pol0), (yy, pol1)):
            append_raw_spectrum(
                out_path, spec,
                utc_time=utc_dt,
                az_deg=float("nan"),
                elev_deg=float("nan"),
                glon_deg=float("nan"),
                glat_deg=float("nan"),
                integration_time_s=float(row["EXPOSURE"]),
                polarization=pol,
                ra_deg=float("nan"),
                dec_deg=float("nan"),
                flag_track=flag_track,
                flag_cal=flag_cal,
                humidity_percent=float("nan"),
                temperature_c=temperature_c,
                pressure_mbar=pressure_mbar,
            )

    hdul.close()
    print(f"  -> {out_path}")
    print(f"     {n_rows} dumps x 2 pol = {n_rows*2} spectra, "
          f"{nchan} channels, beam {beam_number}")
    return out_path


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 fast_fits_to_fro.py <file.fits> [out_dir]")
        sys.exit(1)
    convert_fast_fits(sys.argv[1],
                      sys.argv[2] if len(sys.argv) > 2 else None)