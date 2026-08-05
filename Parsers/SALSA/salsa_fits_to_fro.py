"""
salsa_fits_to_fro.py

Batch converter: SALSA FITS files (single spectrum each) -> one FRO v1.6
session file (multiple spectra).

Workflow:
1. List all .fits files found in ~/Downloads, showing timestamp and
   pointing (GLON/GLAT) read from each header.
2. Let the operator select which files belong to the session to convert
   (printer-style selection syntax, e.g. "1-27" or "3,5,7-12").
3. Verify all selected files share the same frequency axis definition
   (NAXIS1, CRPIX1, CRVAL1, CDELT1) - abort otherwise.
4. Ask interactively for facility / provider / instrument / survey /
   dataset / observation_id / output filename, with sensible defaults
   suggested from the FITS headers. Never assume fixed values that could
   be wrong for a different session.
5. For each selected file (processed in chronological order): read
   header + data, convert data to native byte order, compute Az/El from
   GLON/GLAT + DATE-OBS + Onsala site coordinates (explicitly labeled as
   "computed", never as measured), and append the spectrum to the FRO
   v1.6 session file.
6. Print a final summary via read_fro_summary().
"""

import os
import glob
import warnings
import subprocess
from datetime import datetime, timezone

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.time import Time
import astropy.units as u

from fro_format_v16 import create_fro_v16_structure, append_raw_spectrum, read_fro_summary

# Suppress only the known, harmless FITS padding warning from SALSA files
# (null bytes instead of spaces for header block padding).
warnings.filterwarnings("ignore", message=".*null bytes.*")

DEFAULT_INPUT_DIR = os.path.expanduser("~/FRO/Parsers/SALSA")

# Onsala Space Observatory coordinates, verified via the official Chalmers
# "Travel to Onsala Space Observatory" page (GPS station coordinates).
ONSALA_LOCATION = EarthLocation(
    lat=57.3953 * u.deg,
    lon=11.9255 * u.deg,
    height=10 * u.m,
)


def choose_input_dir():
    """Ask the user for the FITS input directory using a zenity folder dialog.

    Falls back to DEFAULT_INPUT_DIR if zenity is unavailable or the user
    cancels the dialog.
    """
    try:
        result = subprocess.run(
            ["zenity", "--file-selection", "--directory",
             "--title=Select the directory containing the SALSA FITS files",
             f"--filename={DEFAULT_INPUT_DIR}/"],
            capture_output=True, text=True, timeout=300)
    except (FileNotFoundError, subprocess.TimeoutExpired):
        print(f"zenity not available, using default: {DEFAULT_INPUT_DIR}")
        return DEFAULT_INPUT_DIR
    chosen = result.stdout.strip()
    if result.returncode != 0 or not chosen:
        print(f"No directory selected, using default: {DEFAULT_INPUT_DIR}")
        return DEFAULT_INPUT_DIR
    return chosen


def ask(prompt, default):
    """Ask the operator for a value, showing a suggested default."""
    val = input(f"{prompt} [{default}]: ").strip()
    return val if val else default


def parse_selection(selection_str, max_n):
    """Parse a printer-style selection string (e.g. '1-3,5,7-9') into a
    sorted list of unique indices, clipped to the valid range."""
    indices = set()
    for part in selection_str.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_str, end_str = part.split("-", 1)
            start, end = int(start_str), int(end_str)
            indices.update(range(start, end + 1))
        else:
            indices.add(int(part))
    return sorted(i for i in indices if 1 <= i <= max_n)


def read_fits_header_info(filepath):
    """Read the key header fields needed from a SALSA FITS file."""
    with fits.open(filepath) as hdul:
        header = hdul[0].header
        info = {
            "filepath": filepath,
            "telescop": header.get("TELESCOP", "unknown"),
            "origin": header.get("ORIGIN", "unknown"),
            "date_obs": header.get("DATE-OBS", ""),
            "obstime": header.get("OBSTIME", np.nan),
            "naxis1": header.get("NAXIS1"),
            "crpix1": header.get("CRPIX1"),
            "crval1": header.get("CRVAL1"),
            "cdelt1": header.get("CDELT1"),
            "glon": header.get("CRVAL2", np.nan),
            "glat": header.get("CRVAL3", np.nan),
        }
    return info


def build_frequency_axis(info):
    """Reconstruct the frequency axis (Hz) from CRPIX1/CRVAL1/CDELT1."""
    n = int(info["naxis1"])
    crpix1 = info["crpix1"]
    crval1 = info["crval1"]
    cdelt1 = info["cdelt1"]
    pixels = np.arange(1, n + 1)  # FITS pixel numbering is 1-based
    return crval1 + (pixels - crpix1) * cdelt1


def compute_az_el(glon_deg, glat_deg, date_obs_str):
    """Compute Az/El from Galactic coordinates + observation time, using
    the Onsala site location. This is a COMPUTED value, not a measured
    one (SALSA does not provide Az/El in the FITS header)."""
    date_obs_dt = datetime.fromisoformat(date_obs_str)
    obstime = Time(date_obs_dt)
    coord = SkyCoord(l=glon_deg * u.deg, b=glat_deg * u.deg, frame="galactic")
    altaz = coord.transform_to(AltAz(obstime=obstime, location=ONSALA_LOCATION))
    return altaz.az.deg, altaz.alt.deg


def main():
    input_dir = choose_input_dir()
    fits_files = sorted(glob.glob(os.path.join(input_dir, "*.fits")))

    if not fits_files:
        print(f"No FITS files found in {input_dir}.")
        return

    print(f"Found {len(fits_files)} FITS files in {input_dir}:\n")

    file_infos = [read_fits_header_info(fp) for fp in fits_files]

    for idx, info in enumerate(file_infos, start=1):
        date_display = info["date_obs"].replace("T", " ") if info["date_obs"] else "unknown"
        print(f"{idx:3d}. {os.path.basename(info['filepath'])}   "
              f"{date_display}   GLON={info['glon']:.1f}  GLAT={info['glat']:.1f}")

    # Defaults are carried over between correction loop iterations, so the
    # operator only needs to re-type what was actually wrong, not everything.
    default_selection = ""
    default_facility = None
    default_provider = "Onsala Space Observatory"
    default_instrument = None
    default_survey = ""
    default_dataset = ""
    default_observation_id = None
    default_output_filename = None
    observation_id_customized = False
    dataset_customized = False
    output_filename_customized = False


    while True:
        print()
        selection_str = ask(
            "Select files to include in this session (e.g. '1-27' or '3,5,7-12')",
            default_selection
        )
        selected_indices = parse_selection(selection_str, len(file_infos))

        if not selected_indices:
            print("No valid files selected. Please try again.")
            continue

        selected_infos = [file_infos[i - 1] for i in selected_indices]

        # Verify all selected files share the same frequency axis definition.
        ref = selected_infos[0]
        mismatch = False
        for info in selected_infos[1:]:
            if (info["naxis1"] != ref["naxis1"] or
                    info["crpix1"] != ref["crpix1"] or
                    info["crval1"] != ref["crval1"] or
                    info["cdelt1"] != ref["cdelt1"]):
                print("ERROR: selected files do not share the same frequency "
                      "axis definition (NAXIS1/CRPIX1/CRVAL1/CDELT1 differ).")
                print(f"  Reference file: {ref['filepath']}")
                print(f"  Mismatched file: {info['filepath']}")
                print("Please correct your selection.")
                mismatch = True
                break
        if mismatch:
            default_selection = selection_str
            continue

        print(f"\n{len(selected_infos)} files selected, frequency axis consistent.\n")

        # Sort selected files chronologically by DATE-OBS.
        selected_infos.sort(key=lambda i: i["date_obs"])

        if default_instrument is None:
            default_instrument = ref["telescop"].strip().capitalize()
        if default_facility is None:
            default_facility = ref["origin"]
        auto_observation_id = (
            f"SALSA_{selected_infos[0]['date_obs'][:10].replace('-', '')}_session"
        )
        auto_output_filename = (
            f"SALSA_{selected_infos[0]['date_obs'][:10].replace('-', '')}_session.fro"
        )
        glon_values = [i["glon"] for i in selected_infos]
        glat_values = [i["glat"] for i in selected_infos]
        auto_dataset = (
            f"GLon {min(glon_values):.1f} to {max(glon_values):.1f}, "
            f"GLat {min(glat_values):.1f} to {max(glat_values):.1f}"
        )
        if not observation_id_customized:
            default_observation_id = auto_observation_id
        if not output_filename_customized:
            default_output_filename = auto_output_filename
        if not dataset_customized:
            default_dataset = auto_dataset

        facility = ask("Facility", default_facility)
        provider = ask("Provider", default_provider)
        instrument = ask("Instrument", default_instrument)
        survey = ask("Survey", default_survey)
        dataset = ask("Dataset", default_dataset)
        if dataset != auto_dataset:
            dataset_customized = True
        observation_id = ask("Observation ID", default_observation_id)
        output_filename = ask("Output .fro filename", default_output_filename)
        if observation_id != auto_observation_id:
            observation_id_customized = True
        if output_filename != auto_output_filename:
            output_filename_customized = True

        # Remember current values as defaults in case we loop back.
        default_selection = selection_str
        default_facility = facility
        default_provider = provider
        default_instrument = instrument
        default_survey = survey
        default_dataset = dataset
        default_observation_id = observation_id
        default_output_filename = output_filename

        overwrite_ok = True
        if os.path.exists(output_filename):
            while True:
                confirm = ask(
                    f"File '{output_filename}' already exists. Overwrite? (y/N)", "N"
                ).strip().lower()
                if confirm in ("y", "yes"):
                    break
                if confirm in ("n", "no", ""):
                    print("File not overwritten. Please correct the filename or selection.")
                    overwrite_ok = False
                    break
                print(f"'{confirm}' is not a valid answer. Please type y or n.")
        if not overwrite_ok:
            continue

        print("\nPlease review the values entered before writing the file:")
        print(f"  Selected files:  {selection_str} ({len(selected_infos)} files)")
        print(f"  Facility:        {facility}")
        print(f"  Provider:        {provider}")
        print(f"  Instrument:      {instrument}")
        print(f"  Survey:          {survey}")
        print(f"  Dataset:         {dataset}")
        print(f"  Observation ID:  {observation_id}")
        print(f"  Output filename: {output_filename}")

        confirmed = False
        while True:
            proceed = ask("Is everything correct? (y/N)", "N").strip().lower()
            if proceed in ("y", "yes"):
                confirmed = True
                break
            if proceed in ("n", "no", ""):
                print("Let's correct the values.")
                break
            print(f"'{proceed}' is not a valid answer. Please type y or n.")

        if confirmed:
            break  # exit the correction loop, proceed to write the file
        # else: loop back to selection/metadata prompts, keeping current values as defaults

    frequency_axis_hz = build_frequency_axis(ref)
    import_date = datetime.now(timezone.utc).isoformat()

    create_fro_v16_structure(
        output_filename,
        frequency_axis_hz=frequency_axis_hz,
        facility=facility,
        instrument=instrument,
        observation_id=observation_id,
        survey=survey,
        dataset=dataset,
        provider=provider,
        import_date=import_date,
        observatory_name=f"Onsala Space Observatory - {instrument}",
        latitude_deg=ONSALA_LOCATION.lat.deg,
        longitude_deg=ONSALA_LOCATION.lon.deg,
        elevation_m=ONSALA_LOCATION.height.to(u.m).value,
        receiver="SALSA remote receiver",
        antenna=f"{instrument} (2.3m dish)",
        lna="n/a (remote facility, chain not controlled by FRO)",
        center_frequency_hz=float(ref["crval1"]),
        bandwidth_hz=abs(float(ref["naxis1"]) * float(ref["cdelt1"])),
        integration_time_s=np.nan,  # varies per spectrum, see Time/integration_time_s
        notes=f"Converted from {len(selected_infos)} SALSA FITS files via salsa_fits_to_fro.py. Az/El computed from GLON/GLAT + DATE-OBS via Astropy (compute_az_el), not telescope telemetry.",
        processed_by="Francesco Radio Observatory",
        pointing_source="computed",
    )

    print(f"\nCreated {output_filename}, writing {len(selected_infos)} spectra...\n")

    for info in selected_infos:
        with fits.open(info["filepath"]) as hdul:
            data = hdul[0].data

        spectrum = data.reshape(-1).astype(np.float32)  # converts to native byte order

        az_deg, el_deg = compute_az_el(info["glon"], info["glat"], info["date_obs"])
        utc_time = datetime.fromisoformat(info["date_obs"])

        append_raw_spectrum(
            output_filename,
            spectrum,
            utc_time=utc_time,
            az_deg=az_deg,
            elev_deg=el_deg,
            glon_deg=info["glon"],
            glat_deg=info["glat"],
            integration_time_s=info["obstime"],
        )

        print(f"  {os.path.basename(info['filepath'])}: "
              f"GLON={info['glon']:.1f} GLAT={info['glat']:.1f} -> "
              f"Az={az_deg:.1f} El={el_deg:.1f} (integration {info['obstime']:.1f}s)")

    print("\nDone.\n\n" + "="*42 + "\nSummary:\n" + "="*42)
    read_fro_summary(output_filename)
    print("="*42)


if __name__ == "__main__":
    main()