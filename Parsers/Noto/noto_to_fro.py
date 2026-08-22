#!/usr/bin/env python3
# =============================================================================
# noto_to_fro.py
#
# Parser: Noto 32m fitszilla FITS file -> FRO format v17 (HDF5)
#
# Project: FRO (Francesco Radio Observatory) / HALO
# Author: Francesco Di Giovanni (Bolzano, Italy)
# Created: 2026-08-21 (Parte 67)
#
# Reads a fitszilla file (Italian radiotelescopes standard, e.g. Noto 32m)
# and converts it to FRO v17 HDF5 using fro_format_v17.py.
#
# Design rules:
#   - Fully dynamic: every value is read from the FITS file itself
#     (no hardcoded frequencies, bin counts, site coordinates, etc.),
#     because future fitszilla versions may differ (see HISTORY cards
#     in the Primary header).
#   - Nothing from the source file is discarded: full Stokes (LL, RR,
#     cross Re/Im), per-row flags, weather, calibration marks and scan
#     metadata are all preserved (FRO v17 groups Quality, Environment,
#     Calibration).
#   - Frequency axis: fitszilla SECTION TABLE reports the backend
#     sampled band (frequency=0, bandWidth=total). The sky frequency
#     axis is reconstructed from RF INPUTS localOscillator and the
#     sign of RF INPUTS bandWidth:
#         bandWidth < 0  ->  LO above the band, spectrum inverted:
#                            f_sky(bin) = LO - bin * (sect_bw / bins)
#                            (axis and spectra are then flipped so the
#                             stored axis is ascending)
#         bandWidth >= 0 ->  f_sky(bin) = rf_freq + bin * (sect_bw / bins)
#     The negative-bandWidth recipe was verified empirically on the
#     L-band reference file (passband envelope, GSM1800, GNSS, radar,
#     protected band all at the correct frequencies).
#     NOTE: srttools frequency axis is NOT used here (verified to
#     produce an impossible 2300-3800 MHz axis for this file).
# =============================================================================

import sys
import os
import argparse
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u

# fro_format_v17.py must be in the same directory as this script
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fro_format_v17 import (
    create_fro_v17_structure,
    append_raw_spectrum,
    append_cross_spectrum,
)


def rad_to_deg(x):
    """Convert radians to degrees (fitszilla stores angles in radians)."""
    return float(np.degrees(x))


def build_frequency_axis(section_row, rf_rows):
    """
    Build the ascending sky-frequency axis (Hz) for one section.

    Returns (freq_axis_hz, inverted) where 'inverted' is True when the
    raw spectra are stored with descending sky frequency and must be
    flipped to match the ascending axis.
    """
    bins = int(section_row["bins"])
    sect_bw_mhz = float(section_row["bandWidth"])   # backend sampled band
    step_mhz = sect_bw_mhz / bins

    # All RF inputs of a section share the same LO and bandWidth sign;
    # verify instead of assuming.
    lo_values = set(float(r["localOscillator"]) for r in rf_rows)
    bw_values = [float(r["bandWidth"]) for r in rf_rows]
    if len(lo_values) != 1:
        raise ValueError(f"RF inputs of one section have different LOs: {lo_values}")
    lo_mhz = lo_values.pop()

    if all(bw < 0 for bw in bw_values):
        # LO above the band: spectrum inverted (verified recipe)
        f_mhz = lo_mhz - np.arange(bins) * step_mhz
        inverted = True
    elif all(bw >= 0 for bw in bw_values):
        rf_freq_mhz = float(rf_rows[0]["frequency"])
        f_mhz = rf_freq_mhz + np.arange(bins) * step_mhz
        inverted = False
        print("WARNING: non-negative RF bandWidth branch is untested "
              "(reference file has bandWidth < 0). Verify the axis.")
    else:
        raise ValueError("Mixed bandWidth signs among RF inputs of one section")

    if inverted:
        f_mhz = f_mhz[::-1]   # make the stored axis ascending

    return f_mhz * 1e6, inverted


def convert(fits_path, out_path=None):
    """Convert one fitszilla file to FRO v17 HDF5."""
    hdul = fits.open(fits_path)
    primary = hdul[0].header
    section_tbl = hdul["SECTION TABLE"].data
    rf_tbl = hdul["RF INPUTS"].data
    data_tbl = hdul["DATA TABLE"].data

    # ---------------- Observatory (from Primary, radians/meters) ----------
    site_lat_deg = rad_to_deg(primary["SiteLatitude"])
    site_lon_deg = rad_to_deg(primary["SiteLongitude"])
    site_height_m = float(primary["SiteHeight"])
    antenna_name = str(primary.get("ANTENNA", "unknown")).strip()
    receiver_code = str(primary.get("Receiver Code", "unknown")).strip()

    # ---------------- Observation metadata ---------------------------------
    source_name = str(primary.get("SOURCE", "")).strip()
    schedule_name = str(primary.get("ScheduleName", "")).strip()
    scan_id = str(primary.get("SCANID", ""))
    subscan_id = str(primary.get("SubScanID", ""))
    subscan_type = str(primary.get("SubScanType", "")).strip()
    signal = str(primary.get("SIGNAL", "")).strip()
    # fitszilla VLSR unit is assumed to be km/s (srttools/MBFITS convention);
    # FRO v17 stores m/s. Value is 0.0 in the reference file.
    vlsr_mps = float(primary.get("VLSR", np.nan)) * 1000.0

    # ---------------- Section / RF inputs (only one section expected here,
    # but iterate dynamically: data columns are named Ch<id>) ---------------
    if len(section_tbl) != 1:
        print(f"WARNING: {len(section_tbl)} sections found; "
              "this parser currently handles the first one only.")
    section_row = section_tbl[0]
    section_id = int(section_row["id"])
    bins = int(section_row["bins"])
    rf_rows = [r for r in rf_tbl if int(r["section"]) == section_id]
    if len(rf_rows) == 0:
        raise ValueError(f"No RF inputs found for section {section_id}")

    polarizations = [str(r["polarization"]).strip() for r in rf_rows]
    cal_marks_k = [float(r["calibrationMark"]) for r in rf_rows]
    attenuations_db = [float(r["attenuation"]) for r in rf_rows]

    freq_axis_hz, inverted = build_frequency_axis(section_row, rf_rows)
    center_freq_hz = float(np.mean(freq_axis_hz))
    bandwidth_hz = abs(float(section_row["bandWidth"])) * 1e6
    sample_rate_hz = float(section_row["sampleRate"]) * 1e6

    # ---------------- Per-row time axis ------------------------------------
    mjd = np.asarray(data_tbl["time"], dtype=np.float64)
    times_utc = Time(mjd, format="mjd", scale="utc")
    # Integration time is not in the Primary header: derive it from the
    # median spacing of the DATA TABLE time column (days -> seconds).
    if len(mjd) > 1:
        integration_s = float(np.median(np.diff(mjd)) * 86400.0)
    else:
        integration_s = np.nan

    # ---------------- Output file ------------------------------------------
    if out_path is None:
        base = os.path.splitext(os.path.basename(fits_path))[0]
        out_path = os.path.join(os.path.dirname(os.path.abspath(fits_path)),
                                base + "_FRO_v17.h5")

    observation_id = (f"{antenna_name}_{os.path.splitext(schedule_name)[0]}"
                      f"_scan{scan_id}_sub{subscan_id}")

    create_fro_v17_structure(
        out_path,
        frequency_axis_hz=freq_axis_hz,
        facility=antenna_name,
        instrument=receiver_code,
        observation_id=observation_id,
        survey="",
        dataset="",
        provider="private test file (not for distribution)",
        import_date=datetime.now(timezone.utc).isoformat(),
        observatory_name=antenna_name,
        latitude_deg=site_lat_deg,
        longitude_deg=site_lon_deg,
        elevation_m=site_height_m,
        receiver=receiver_code,
        antenna=antenna_name,
        lna="n/a",
        center_frequency_hz=center_freq_hz,
        sample_rate_hz=sample_rate_hz,
        bandwidth_hz=bandwidth_hz,
        integration_time_s=integration_s,
        notes=(f"Converted from fitszilla file '{os.path.basename(fits_path)}'. "
               f"Spectrum inverted in source file: {inverted}. "
               "Frequency axis rebuilt from RF INPUTS localOscillator and "
               "SECTION bandWidth/bins (empirically verified recipe)."),
        temperature_calibration_note=(
            "Raw backend counts (INT32 -> float32); no temperature "
            "calibration applied. Calibration mark temperatures per RF "
            "input are stored in the Calibration group."),
        processed_by="noto_to_fro.py",
        pointing_source="measured",
        target=source_name,
        vlsr_mps=vlsr_mps,
        scan_id=scan_id,
        subscan_id=subscan_id,
        subscan_type=subscan_type,
        signal=signal,
        schedule_name=schedule_name,
        cal_polarizations=polarizations,
        cal_mark_temp_k=cal_marks_k,
        attenuation_db=attenuations_db,
        include_cross=True,
    )

    # ---------------- Per-row conversion ------------------------------------
    ch_col = f"Ch{section_id}"
    n_rows = len(data_tbl)
    for i in range(n_rows):
        row = data_tbl[i]
        # Stokes layout: 4 consecutive blocks of <bins> INT32 (LL, RR, Re, Im)
        stokes = np.asarray(row[ch_col], dtype=np.int32).reshape(4, bins)
        ll = stokes[0].astype(np.float32)
        rr = stokes[1].astype(np.float32)
        cross_re = stokes[2].astype(np.float32)
        cross_im = stokes[3].astype(np.float32)
        if inverted:
            ll = ll[::-1]
            rr = rr[::-1]
            cross_re = cross_re[::-1]
            cross_im = cross_im[::-1]

        az_deg = rad_to_deg(row["az"])
        el_deg = rad_to_deg(row["el"])
        ra_deg = rad_to_deg(row["raj2000"])
        dec_deg = rad_to_deg(row["decj2000"])
        # Galactic coordinates computed from measured RA/Dec (J2000)
        gal = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg,
                       frame="icrs").galactic
        glon_deg = float(gal.l.deg)
        glat_deg = float(gal.b.deg)

        utc_iso = times_utc[i].to_datetime(timezone=timezone.utc)
        flag_cal = int(row["flag_cal"])
        flag_track = int(row["flag_track"])
        humidity, temperature, pressure = [float(x) for x in row["weather"]]

        # Two rows per integration (one per polarization), duplicated
        # metadata on both (intentional redundancy, decided in Parte 30)
        for spec, pol in ((ll, polarizations[0]), (rr, polarizations[1])):
            append_raw_spectrum(
                out_path, spec,
                utc_time=utc_iso,
                az_deg=az_deg, elev_deg=el_deg,
                glon_deg=glon_deg, glat_deg=glat_deg,
                integration_time_s=integration_s,
                polarization=pol,
                ra_deg=ra_deg, dec_deg=dec_deg,
                flag_track=flag_track, flag_cal=flag_cal,
                humidity_percent=humidity,
                temperature_c=temperature,
                pressure_mbar=pressure,
            )
        # One cross-spectrum row per integration
        append_cross_spectrum(out_path, cross_re, cross_im)

    hdul.close()

    # ---------------- Summary ------------------------------------------------
    print("=" * 60)
    print("Noto fitszilla -> FRO v17 conversion complete")
    print("=" * 60)
    print(f"Input file      : {fits_path}")
    print(f"Output file     : {out_path}")
    print(f"Facility        : {antenna_name}  (receiver {receiver_code})")
    print(f"Source          : {source_name}")
    print(f"Schedule        : {schedule_name}")
    print(f"Scan/Subscan    : {scan_id}/{subscan_id}  type={subscan_type}  signal={signal}")
    print(f"Integrations    : {n_rows}  x  {integration_s*1000:.1f} ms")
    print(f"Rows written    : {n_rows * len(polarizations)} spectra "
          f"({' + '.join(polarizations)}) + {n_rows} cross rows")
    print(f"Channels        : {bins}")
    print(f"Frequency axis  : {freq_axis_hz[0]/1e6:.3f} - {freq_axis_hz[-1]/1e6:.3f} MHz "
          f"(step {abs(freq_axis_hz[1]-freq_axis_hz[0])/1e3:.3f} kHz, "
          f"inverted in source: {inverted})")
    print(f"Site            : lat {site_lat_deg:.4f}  lon {site_lon_deg:.4f}  "
          f"alt {site_height_m:.0f} m")
    print("=" * 60)
    return out_path


def main():
    parser = argparse.ArgumentParser(
        description="Convert a Noto 32m fitszilla FITS file to FRO format v17 (HDF5)")
    parser.add_argument("fits_file", help="input fitszilla .fits file")
    parser.add_argument("-o", "--output", default=None,
                        help="output .h5 path (default: alongside input)")
    args = parser.parse_args()
    convert(args.fits_file, args.output)


if __name__ == "__main__":
    main()