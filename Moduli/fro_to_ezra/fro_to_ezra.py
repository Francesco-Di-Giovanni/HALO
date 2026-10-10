"""
fro_to_ezra.py

Converts a FRO file (v1.2 or v1.3, HDF5) into the plain-text format
accepted by ezRA's ezCon.py.

Usage:
    python3 fro_to_ezra.py <input.fro> [output.txt]

If output.txt is omitted, it defaults to <input>_ezra.txt in the same
directory as the input file.

Key design notes:
- Reads facility name and observatory coordinates dynamically from the
  /Observatory group (no hardcoded values).
- Handles multiple spectra per file (Spectra/raw_spectrum is a 2D array:
  N_spectra x channels), not just a single spectrum.
- Uses the real per-spectrum UTC timestamps stored in Time/utc, instead
  of generating fake ones.
- Reads per-spectrum azimuth/elevation from /Pointing if present
  (FRO v1.3+). If the group is missing (older v1.2 files), az/el are
  treated as NaN for all spectra. A new "az ... el ..." pointing line
  is written to the output only when the value changes from the
  previous spectrum, matching ezCon.py's expected format (it keeps
  using the last pointing line until a new one appears).
- All spectra are currently flagged "a" (antenna/on-source). No
  REF/calibration distinction is made yet, since no hot/cold or
  noise-diode calibration procedure is implemented at the acquisition
  stage as of this version.
"""

import sys
import os
import h5py
import numpy as np
from datetime import datetime


def parse_fro_timestamp(iso_string):
    """
    Parses an ISO-format timestamp string (as stored in Time/utc,
    possibly with a UTC offset) into a datetime object.
    """
    return datetime.fromisoformat(iso_string)


def read_fro_file(input_fro):
    """
    Reads all data needed from a FRO file.
    Returns a dict with observatory info, frequency axis, spectra,
    timestamps, and pointing (azimuth/elevation) arrays.
    """
    with h5py.File(input_fro, "r") as f:
        observatory_name = f["Observatory"].attrs.get("name", "Unknown")
        lat = float(f["Observatory"].attrs.get("latitude_deg", np.nan))
        lon = float(f["Observatory"].attrs.get("longitude_deg", np.nan))
        amsl = float(f["Observatory"].attrs.get("elevation_m", np.nan))

        freq_hz = f["Spectra/frequency_axis_hz"][:]
        spectra = f["Spectra/raw_spectrum"][:]  # shape: (N_spectra, channels)

        n_spectra = spectra.shape[0]

        utc_raw = f["Time/utc"][:]
        timestamps = [parse_fro_timestamp(t.decode() if isinstance(t, bytes) else t)
                      for t in utc_raw]

        if "Pointing" in f:
            az = f["Pointing/azimuth_deg"][:]
            elev = f["Pointing/elevation_deg"][:]
        else:
            az = np.full(n_spectra, np.nan)
            elev = np.full(n_spectra, np.nan)

    return {
        "observatory_name": observatory_name,
        "lat": lat,
        "lon": lon,
        "amsl": amsl,
        "freq_hz": freq_hz,
        "spectra": spectra,
        "timestamps": timestamps,
        "az": az,
        "elev": elev,
    }


def write_ezra_file(data, output_txt):
    """
    Writes the ezRA-format text file from the data read from a FRO file.
    """
    freq_mhz = data["freq_hz"] / 1e6
    freq_min = freq_mhz[0]
    freq_max = freq_mhz[-1]
    freq_bin_qty = len(freq_mhz)
    n_spectra = data["spectra"].shape[0]

    with open(output_txt, "w") as out:
        # REQUIRED: ezCon only accepts files whose first line starts with
        # "from ezCol" (ezCon.py, check fileLine[:10] != 'from ezCol').
        # Do NOT rename to fro_to_ezra.py or ezCon will silently skip the file.
        out.write("from ezColFRO.py\n")
        out.write(
            f"lat {data['lat']} long {data['lon']} amsl {data['amsl']} "
            f"name {data['observatory_name']} ezb\n"
        )
        out.write(f"freqMin {freq_min:.6f} freqMax {freq_max:.6f} freqBinQty {freq_bin_qty}\n")
        out.write("# times are in UTC\n")

        last_az = None
        last_elev = None

        for i in range(n_spectra):
            az_i = data["az"][i]
            elev_i = data["elev"][i]

            # Write a new pointing line only if it changed since the
            # previous spectrum (or this is the first spectrum).
            az_changed = last_az is None or not (
                (np.isnan(az_i) and np.isnan(last_az)) or az_i == last_az
            )
            elev_changed = last_elev is None or not (
                (np.isnan(elev_i) and np.isnan(last_elev)) or elev_i == last_elev
            )
            if az_changed or elev_changed:
                out.write(f"az {az_i:g} el {elev_i:g}\n")
                last_az = az_i
                last_elev = elev_i

            row_time = data["timestamps"][i].strftime("%Y-%m-%dT%H:%M:%S")
            out.write(row_time)
            for value in data["spectra"][i]:
                out.write(f" {value:.8g}")
            out.write(" a\n")  # all spectra flagged as antenna/on-source for now

    print(f"Written: {output_txt}")
    print(f"Spectra: {n_spectra}")
    print(f"Bins: {freq_bin_qty}")
    print(f"Frequency range MHz: {freq_min:.6f} - {freq_max:.6f}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 fro_to_ezra.py <input.fro> [output.txt]")
        sys.exit(1)

    input_fro = sys.argv[1]

    if len(sys.argv) >= 3:
        output_txt = sys.argv[2]
    else:
        base = os.path.splitext(input_fro)[0]
        output_txt = base + "_ezra.txt"

    data = read_fro_file(input_fro)
    write_ezra_file(data, output_txt)
