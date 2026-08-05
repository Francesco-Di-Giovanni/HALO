#!/usr/bin/env python3
# gbt_sdfits_to_fro.py
#
# Parser: GBT SDFITS -> FRO format v17.
#
# Reference file: TREG_091209.cal.acs.fits (Bob Garwood's reference example
# from the NASA FITS registry, GBTIDL calibrated output "keep").
# Structure verified empirically on 2026-07-13 (Parte 34):
#   * One or more BinTable extensions named "SINGLE DISH", 69 columns.
#   * Each FITS row = ONE spectrum of ONE polarization (PLNUM), one IF
#     window (IFNUM), one feed (FDNUM), one integration (INT) of one SCAN.
#   * Frequency axis is per-row WCS: CRVAL1/CRPIX1/CDELT1 (linear, honest).
#     CDELT1 may be negative -> spectrum and axis are flipped to ascending.
#   * A single extension may contain several IF windows with different
#     center frequencies (verified in Ext 2: IFNUM 0/1). Since FRO v17
#     stores ONE frequency axis per file, output is split:
#     one .fro file per (extension, IFNUM) group.
#
# Environment units verified empirically on this file:
#   TAMBIENT [K], PRESSURE [mmHg], HUMIDITY [fraction 0-1]
# Converted to v17 native units: temperature_c, pressure_mbar (hPa),
# humidity_percent. Conversions are exact and reversible (documented
# in the Notes field).
#
# Honesty notes:
#   * AZIMUTH/ELEVATIO are measured values -> pointing_source="measured".
#   * Galactic coordinates are computed from CRVAL2/CRVAL3 (RA/DEC J2000)
#     via Astropy -> stated in Notes.
#   * SDFITS has no tracking flag -> flag_track = -1 (unknown).
#   * flag_cal = CAL column (True on cal-on phases of raw files).

import os
import sys

import numpy as np
from datetime import datetime, timezone

from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fro_format_v17 import (
    create_fro_v17_structure,
    append_raw_spectrum,
)

# FITS polarization codes (standard convention, stored in CRVAL4)
POL_CODES = {
    1: "I", 2: "Q", 3: "U", 4: "V",
    -1: "RR", -2: "LL", -3: "RL", -4: "LR",
    -5: "XX", -6: "YY", -7: "XY", -8: "YX",
}

MMHG_TO_HPA = 1.3332239


def pol_label(crval4):
    """Map FITS CRVAL4 polarization code to a human-readable label."""
    return POL_CODES.get(int(crval4), f"POL{int(crval4)}")


def freq_axis_from_row(row, n_channels):
    """Rebuild the sky frequency axis [Hz] from per-row WCS keywords.

    FITS convention: pixel index is 1-based, CRPIX1 marks the reference
    pixel where frequency equals CRVAL1.
    """
    crval1 = float(row["CRVAL1"])
    crpix1 = float(row["CRPIX1"])
    cdelt1 = float(row["CDELT1"])
    idx = np.arange(1, n_channels + 1, dtype=np.float64)
    return crval1 + (idx - crpix1) * cdelt1


def convert_extension_group(hdu, ext_number, rows_idx, ifnum, fits_path,
                            out_dir):
    """Convert one (extension, IFNUM) group of rows to a single .fro file."""
    data = hdu.data
    first = data[rows_idx[0]]
    n_channels = len(first["DATA"])

    # ---------------- Frequency axis (from first row of the group) --------
    freq_axis_hz = freq_axis_from_row(first, n_channels)
    cdelt1 = float(first["CDELT1"])
    flipped = cdelt1 < 0.0
    if flipped:
        freq_axis_hz = freq_axis_hz[::-1]

    # Check frequency-axis stability across the group (doppler tracking
    # may shift CRVAL1 slightly row to row). We keep the first row's axis
    # and record the maximum deviation.
    crval1_all = np.array([float(data[i]["CRVAL1"]) for i in rows_idx])
    max_crval1_dev_hz = float(np.max(np.abs(crval1_all - crval1_all[0])))

    # ---------------- Global metadata (from first row) --------------------
    telescope = str(first["TELESCOP"]).strip()
    frontend = str(first["FRONTEND"]).strip()
    backend = str(first["BACKEND"]).strip()
    target = str(first["OBJECT"]).strip()
    projid = str(first["PROJID"]).strip()
    observer = str(first["OBSERVER"]).strip()
    obsmode = str(first["OBSMODE"]).strip()
    scan_id = int(first["SCAN"])
    site_lon_deg = float(first["SITELONG"])
    site_lat_deg = float(first["SITELAT"])
    site_elev_m = float(first["SITEELEV"])
    bandwidth_hz = float(first["BANDWID"])
    center_freq_hz = float(first["OBSFREQ"])
    restfreq_hz = float(first["RESTFREQ"])
    veldef = str(first["VELDEF"]).strip()
    velocity_mps = float(first["VELOCITY"])
    tunit = str(first["TUNIT7"]).strip()
    exposure_s = float(first["EXPOSURE"])

    # Calibration info: one TCAL per polarization present in the group
    pols_in_group = []
    tcal_in_group = []
    for i in rows_idx:
        lbl = pol_label(data[i]["CRVAL4"])
        if lbl not in pols_in_group:
            pols_in_group.append(lbl)
            tcal_in_group.append(float(data[i]["TCAL"]))

    # ---------------- Output file -----------------------------------------
    base = os.path.splitext(os.path.basename(fits_path))[0]
    # Double extension like ".cal.acs.fits" -> strip once more is not done:
    # we keep the base as-is for traceability.
    out_name = f"{base}_ext{ext_number}_if{ifnum}_FRO_v17.h5"
    out_path = os.path.join(out_dir, out_name)

    observation_id = (f"{telescope}_{projid}_scan{scan_id}"
                      f"_ext{ext_number}_if{ifnum}")

    notes = (
        f"Converted from GBT SDFITS file '{os.path.basename(fits_path)}', "
        f"extension {ext_number}, IFNUM {ifnum} "
        f"({len(rows_idx)} rows, {n_channels} channels). "
        f"GBTIDL calibrated output (TUNIT '{tunit}'), not raw sdfits "
        "filler data. "
        f"OBSMODE '{obsmode}'. Observer: {observer}. "
        "Frequency axis rebuilt from per-row WCS (CRVAL1/CRPIX1/CDELT1); "
        f"axis taken from first row, max CRVAL1 deviation across rows: "
        f"{max_crval1_dev_hz:.3f} Hz. "
        f"Spectrum and axis flipped to ascending: {flipped}. "
        f"CTYPE1 '{first['CTYPE1'].strip()}' (topocentric observed "
        f"frequencies); VELDEF '{veldef}', RESTFREQ {restfreq_hz:.1f} Hz. "
        "Galactic coordinates computed from CRVAL2/CRVAL3 (RA/DEC J2000) "
        "via Astropy. AZIMUTH/ELEVATIO are measured values. "
        "Environment converted from SDFITS native units: TAMBIENT [K] -> "
        "temperature_c, PRESSURE [mmHg] -> pressure_mbar (hPa), HUMIDITY "
        "[fraction 0-1] -> humidity_percent (exact, reversible). "
        "flag_track = -1: SDFITS provides no tracking flag (unknown). "
        "flag_cal = CAL column."
    )

    temperature_calibration_note = (
        f"GBTIDL calibrated antenna temperature (TUNIT '{tunit}'). "
        "Per-row TSYS available in source file; first-row TSYS: "
        f"{float(first['TSYS']):.3f} K. TCAL per polarization stored in "
        "the Calibration group."
    )

    create_fro_v17_structure(
        out_path,
        frequency_axis_hz=freq_axis_hz,
        facility=telescope,
        instrument=f"{frontend}/{backend}",
        observation_id=observation_id,
        survey="",
        dataset="",
        provider="NASA FITS registry reference file (Bob Garwood)",
        import_date=datetime.now(timezone.utc).isoformat(),
        observatory_name=telescope,
        latitude_deg=site_lat_deg,
        longitude_deg=site_lon_deg,
        elevation_m=site_elev_m,
        receiver=frontend,
        antenna=telescope,
        lna="n/a",
        center_frequency_hz=center_freq_hz,
        sample_rate_hz=float("nan"),
        bandwidth_hz=bandwidth_hz,
        integration_time_s=exposure_s,
        notes=notes,
        temperature_calibration_note=temperature_calibration_note,
        processed_by="gbt_sdfits_to_fro.py",
        pointing_source="measured",
        target=target,
        vlsr_mps=velocity_mps,
        scan_id=scan_id,
        subscan_id=ifnum,
        subscan_type=obsmode,
        signal="unknown",
        schedule_name=projid,
        cal_polarizations=pols_in_group,
        cal_mark_temp_k=tcal_in_group,
        attenuation_db=[float("nan")] * len(pols_in_group),
        include_cross=False,
    )

    # ---------------- Per-row conversion -----------------------------------
    for i in rows_idx:
        row = data[i]
        spec = np.asarray(row["DATA"], dtype=np.float32)
        if flipped:
            spec = spec[::-1]

        az_deg = float(row["AZIMUTH"])
        el_deg = float(row["ELEVATIO"])
        ra_deg = float(row["CRVAL2"])
        dec_deg = float(row["CRVAL3"])
        gal = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg,
                       frame="icrs").galactic
        glon_deg = float(gal.l.deg)
        glat_deg = float(gal.b.deg)

        # DATE-OBS like '2009-12-10T02:46:39.00' -> datetime (UTC)
        utc_dt = datetime.fromisoformat(
            str(row["DATE-OBS"]).strip()).replace(tzinfo=timezone.utc)

        flag_cal = int(bool(row["CAL"]))
        flag_track = -1  # unknown: SDFITS has no tracking flag

        temperature_c = float(row["TAMBIENT"]) - 273.15
        pressure_hpa = float(row["PRESSURE"]) * MMHG_TO_HPA
        humidity_percent = float(row["HUMIDITY"]) * 100.0

        append_raw_spectrum(
            out_path, spec,
            utc_time=utc_dt,
            az_deg=az_deg, elev_deg=el_deg,
            glon_deg=glon_deg, glat_deg=glat_deg,
            integration_time_s=float(row["EXPOSURE"]),
            polarization=pol_label(row["CRVAL4"]),
            ra_deg=ra_deg, dec_deg=dec_deg,
            flag_track=flag_track, flag_cal=flag_cal,
            humidity_percent=humidity_percent,
            temperature_c=temperature_c,
            pressure_mbar=pressure_hpa,
        )

    print(f"  -> {out_path}")
    print(f"     {len(rows_idx)} spectra, {n_channels} channels, "
          f"target '{target}', pols {pols_in_group}")
    return out_path


def convert_gbt_sdfits(fits_path, out_dir=None):
    """Convert every (SINGLE DISH extension, IFNUM) group to a .fro file."""
    fits_path = os.path.abspath(fits_path)
    if out_dir is None:
        out_dir = os.path.dirname(fits_path)

    hdul = fits.open(fits_path)
    written = []
    for ext_number, hdu in enumerate(hdul):
        if not isinstance(hdu, fits.BinTableHDU):
            continue
        extname = str(hdu.header.get("EXTNAME", "")).strip()
        if extname.upper() != "SINGLE DISH":
            continue

        ifnums = np.asarray(hdu.data["IFNUM"], dtype=int)
        for ifnum in sorted(set(ifnums.tolist())):
            rows_idx = [k for k in range(len(ifnums)) if ifnums[k] == ifnum]
            print(f"Extension {ext_number}, IFNUM {ifnum}: "
                  f"{len(rows_idx)} rows")
            written.append(
                convert_extension_group(hdu, ext_number, rows_idx, ifnum,
                                        fits_path, out_dir))
    hdul.close()

    if not written:
        print("No 'SINGLE DISH' binary table extensions found.")
    return written


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 gbt_sdfits_to_fro.py <file.fits> [out_dir]")
        sys.exit(1)
    convert_gbt_sdfits(sys.argv[1],
                       sys.argv[2] if len(sys.argv) > 2 else None)