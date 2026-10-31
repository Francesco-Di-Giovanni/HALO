"""
salsa_fits_to_fro.py - Parser: SALSA FITS spectra (one per file), batched
into one session -> FRO v1.8 (fro_format_v1_8).

Every value the FITS header provides is kept; values are computed only when
the header does not provide them, and that is declared in the .fro.

FITS keywords used (checked on SALSA files of July and October 2026):
  DATE-OBS, OBSTIME            time and integration time
  NAXIS1, CRPIX1, CRVAL1, CDELT1  frequency axis (CRVAL1 at pixel CRPIX1)
  CRVAL2, CRVAL3               galactic longitude / latitude
  AZIMUTH, ELEVATIO            commanded Az/El at the start of the
                               integration (missing in files before ~10 July
                               2026: then computed with Astropy, declared)
  VELOSYS, SPECSYS             LSR velocity correction (m/s) and its frame
  RESTFRQ                      always the HI rest frequency, whatever the
                               observed line: kept as source_restfrq_hz, the
                               observed line is chosen at import
  ORIGIN, TELESCOP             facility and instrument
  (VELO-LSR is the same value as VELOSYS in km/s, legacy: not stored twice)

Workflow:
1. Folder with the SALSA FITS files: command-line argument, or a zenity
   folder dialog (default ~/FRO/Parsers/SALSA).
2. List the files with timestamp, pointing, VELOSYS and whether Az/El are in
   the header.
3. Questions, one at a time (type < and Enter to go back one step): files
   to include (printer-style, e.g. "1-27" or "3,5,7-12"; same frequency
   axis), observed line (numbered menu), target, survey, output filename.
   Facility, provider and instrument come from the FITS.
4. Overwrite confirmation if the file exists; numbered review: answer n and
   choose the number of the value to correct.
5. Write: one row per FITS file in chronological order.
6. Summary via read_fro_summary().
"""

import os
import glob
import json
import warnings
import subprocess
import sys
from datetime import datetime, timezone

import numpy as np
import h5py
sys.path.insert(0, os.path.expanduser("~/FRO/Moduli/fro_format"))
from astropy.io import fits
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, FK4
from astropy.time import Time
import astropy.units as u

from fro_format_v1_8 import create_fro_structure, append_raw_spectrum, read_fro_summary


KNOWN_TARGETS = os.path.expanduser("~/FRO/Moduli/fro_onoff_stack/known_targets.json")
TARGET_TOL_DEG = 0.5


def proposed_target_survey(sel):
    """Proposal for target and survey: first from a .fro already in the
    folder of the FITS files, then from known_targets.json by coordinates."""
    folder = os.path.dirname(sel[0]["filepath"])
    for fn in sorted(glob.glob(os.path.join(folder, "*.fro"))):
        found = {}

        def grab(obj):
            for k in ("target", "survey"):
                v = obj.attrs.get(k)
                if isinstance(v, bytes):
                    v = v.decode()
                if v is not None and str(v) and k not in found:
                    found[k] = str(v)

        try:
            with h5py.File(fn, "r") as h:
                grab(h)
                h.visititems(lambda n, o: grab(o))
        except OSError:
            continue
        if found:
            return (found.get("target", ""), found.get("survey", ""),
                    os.path.basename(fn))
    try:
        with open(KNOWN_TARGETS) as fh:
            known = json.load(fh)
    except (OSError, ValueError):
        return "", "", ""
    for name, t in known.items():
        if "glon_deg" not in t or "glat_deg" not in t:
            continue
        for i in sel:
            dl = (i["glon"] - t["glon_deg"] + 180.0) % 360.0 - 180.0
            if np.hypot(dl, i["glat"] - t["glat_deg"]) <= TARGET_TOL_DEG:
                return name, t.get("survey", ""), "known_targets.json"
    return "", "", ""

# Suppress only the known, harmless FITS padding warning from SALSA files
# (null bytes instead of spaces for header block padding).
warnings.filterwarnings("ignore", message=".*null bytes.*")

DEFAULT_INPUT_DIR = os.path.expanduser("~/FRO/Parsers/SALSA")
DEFAULT_PROVIDER = "Onsala Space Observatory"

# Onsala Space Observatory coordinates, verified via the official Chalmers
# "Travel to Onsala Space Observatory" page (GPS station coordinates).
# Used only to compute Az/El for files without AZIMUTH/ELEVATIO.
ONSALA_LOCATION = EarthLocation(
    lat=57.3953 * u.deg,
    lon=11.9255 * u.deg,
    height=10 * u.m,
)

# Rest frequencies (Hz) offered in the observed-line menu.
LINES = [
    ("HI 21 cm", 1420.405752e6),
    ("OH 1612 (satellite F=1-2)", 1612.2310e6),
    ("OH 1665 (main F=1-1)", 1665.4018e6),
    ("OH 1667 (main F=2-2)", 1667.3590e6),
    ("OH 1720 (satellite F=2-1)", 1720.5300e6),
]

DOPPLER_COMPUTED_NOTE = (
    "vlsr_correction_computed_mps: LSRK velocity correction computed by "
    "salsa_fits_to_fro.py, independently of the FITS: Astropy barycentric "
    "correction for Onsala (57.3953 N, 11.9255 E, 10 m) at DATE-OBS (start of "
    "integration), plus the solar motion of 20 km/s towards the standard LSRK "
    "apex RA 18h, Dec +30 deg (B1900). Same sign convention as VELOSYS.")

BACK = "<"


def choose_input_dir():
    """Folder from the command line, else a zenity folder dialog, else the
    default folder."""
    if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
        return os.path.abspath(sys.argv[1])
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
    """Ask for a value showing the default; returns BACK for '<'."""
    val = input(f"{prompt} [{default}]: ").strip()
    if val == BACK:
        return BACK
    return val if val else default


def ask_yes_no(prompt):
    ans = ""
    while ans not in ("y", "yes", "n", "no"):
        ans = input(f"{prompt} (y/N) [N]: ").strip().lower() or "n"
        if ans not in ("y", "yes", "n", "no"):
            print(f"  '{ans}' is not a valid answer. Please type y or n.")
    return ans in ("y", "yes")


def parse_selection(selection_str, max_n):
    """Parse a printer-style selection string (e.g. '1-3,5,7-9') into a
    sorted list of unique indices, clipped to the valid range. Raises
    ValueError on malformed input."""
    indices = set()
    for part in selection_str.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start_str, end_str = part.split("-", 1)
            indices.update(range(int(start_str), int(end_str) + 1))
        else:
            indices.add(int(part))
    return sorted(i for i in indices if 1 <= i <= max_n)


def read_fits_header_info(filepath):
    """Read the header fields used from a SALSA FITS file."""
    with fits.open(filepath) as hdul:
        h = hdul[0].header
        info = {
            "filepath": filepath,
            "telescop": h.get("TELESCOP", "unknown"),
            "origin": h.get("ORIGIN", "unknown"),
            "date_obs": h.get("DATE-OBS", ""),
            "obstime": h.get("OBSTIME", np.nan),
            "naxis1": h.get("NAXIS1"),
            "crpix1": h.get("CRPIX1"),
            "crval1": h.get("CRVAL1"),
            "cdelt1": h.get("CDELT1"),
            "glon": h.get("CRVAL2", np.nan),
            "glat": h.get("CRVAL3", np.nan),
            "az": float(h.get("AZIMUTH", np.nan)),
            "el": float(h.get("ELEVATIO", np.nan)),
            "velosys": float(h.get("VELOSYS", np.nan)),
            "specsys": str(h.get("SPECSYS", "")).strip(),
            "restfrq": float(h.get("RESTFRQ", np.nan)),
        }
    info["has_azel"] = bool(np.isfinite(info["az"]) and np.isfinite(info["el"]))
    return info


def build_frequency_axis(info):
    """Reconstruct the frequency axis (Hz) from CRPIX1/CRVAL1/CDELT1."""
    n = int(info["naxis1"])
    pixels = np.arange(1, n + 1)  # FITS pixel numbering is 1-based
    return info["crval1"] + (pixels - info["crpix1"]) * info["cdelt1"]


def compute_az_el(glon_deg, glat_deg, date_obs_str):
    """Az/El from Galactic coordinates + observation time at Onsala. A
    COMPUTED value, used only for files without AZIMUTH/ELEVATIO."""
    obstime = Time(datetime.fromisoformat(date_obs_str))
    coord = SkyCoord(l=glon_deg * u.deg, b=glat_deg * u.deg, frame="galactic")
    altaz = coord.transform_to(AltAz(obstime=obstime, location=ONSALA_LOCATION))
    return altaz.az.deg, altaz.alt.deg


_LSRK_APEX = None


def compute_vlsr_correction(glon_deg, glat_deg, date_obs_str):
    """LSRK velocity correction (m/s) at the start of the integration,
    computed independently of the FITS (see DOPPLER_COMPUTED_NOTE)."""
    global _LSRK_APEX
    if _LSRK_APEX is None:
        apex = SkyCoord(ra=270 * u.deg, dec=30 * u.deg,
                        frame=FK4(equinox="B1900")).icrs
        _LSRK_APEX = 20.0 * apex.cartesian.xyz.value          # km/s
    c = SkyCoord(l=glon_deg * u.deg, b=glat_deg * u.deg, frame="galactic")
    t = Time(datetime.fromisoformat(date_obs_str))
    vb = c.radial_velocity_correction(kind="barycentric", obstime=t,
                                      location=ONSALA_LOCATION).to(u.km / u.s).value
    return float((vb + np.dot(_LSRK_APEX, c.icrs.cartesian.xyz.value)) * 1000.0)


# ---------------------------------------------------------------------------
# Interactive steps. Each step reads/writes the shared state dict and
# returns BACK to go to the previous step, None to go on.
# ---------------------------------------------------------------------------

def make_steps(file_infos):

    def step_selection(st):
        while True:
            s = ask("Select files to include in this session "
                    "(e.g. '1-27' or '3,5,7-12')", st["selection"])
            if s == BACK:
                return BACK
            try:
                idx = parse_selection(s, len(file_infos))
            except ValueError:
                print("  Not a valid selection, try again.")
                continue
            if not idx:
                print("  No valid files selected, try again.")
                continue
            sel = [file_infos[i - 1] for i in idx]
            ref = sel[0]
            key = lambda i: (i["naxis1"], i["crpix1"], i["crval1"], i["cdelt1"])
            bad = [i for i in sel[1:] if key(i) != key(ref)]
            if bad:
                print("  ERROR: the selected files do not share the same "
                      "frequency axis (NAXIS1/CRPIX1/CRVAL1/CDELT1 differ).")
                print(f"  Reference file: {os.path.basename(ref['filepath'])}")
                print(f"  Mismatched file: {os.path.basename(bad[0]['filepath'])}")
                continue
            sel.sort(key=lambda i: i["date_obs"])
            st["selection"] = s
            st["selected"] = sel
            st["axis"] = build_frequency_axis(ref)
            if not st.get("target") and not st.get("survey"):
                t, sv, src = proposed_target_survey(sel)
                if t or sv:
                    st["target"], st["survey"] = t, sv
                    print(f"  Proposed from {src}: target '{t}', survey '{sv}' "
                          "(you can change them in the next steps).")
            if not st["custom"].get("instrument"):
                st["instrument"] = ref["telescop"].strip().capitalize()
            if not st["custom"].get("facility"):
                st["facility"] = ref["origin"].strip()
            n_azel = sum(i["has_azel"] for i in sel)
            print(f"  {len(sel)} files selected, frequency axis consistent; "
                  f"Az/El in the header: {n_azel} of {len(sel)}.")
            return None

    def step_line(st):
        axis = st["axis"]
        fmin, fmax = min(axis[0], axis[-1]), max(axis[0], axis[-1])
        inband = [k for k, (_, f) in enumerate(LINES) if fmin <= f <= fmax]
        n_other = len(LINES) + 1
        print(f"\nBand of the selected files: {fmin / 1e6:.3f} - "
              f"{fmax / 1e6:.3f} MHz")
        rf = st["selected"][0]["restfrq"]
        if np.isfinite(rf):
            print(f"RESTFRQ written in the FITS: {rf / 1e6:.6f} MHz "
                  "(SALSA writes HI whatever the observed line)")
        print("Observed line:")
        for k, (name, f) in enumerate(LINES, 1):
            mark = "   <- in band" if (k - 1) in inband else ""
            print(f"  {k}. {name:28s} {f / 1e6:12.6f} MHz{mark}")
        print(f"  {n_other}. other (type the rest frequency)")
        if st["line_choice"] is None:
            st["line_choice"] = str(inband[0] + 1) if inband else "1"
        while True:
            v = ask("Choose a number", st["line_choice"])
            if v == BACK:
                return BACK
            if v.isdigit() and 1 <= int(v) <= len(LINES):
                st["line_choice"] = v
                st["line_name"], st["rest_hz"] = LINES[int(v) - 1]
                if not fmin <= st["rest_hz"] <= fmax:
                    print("  Note: this line is outside the band of the "
                          "selected files.")
                return None
            if v == str(n_other):
                prev = (f"{st['rest_hz'] / 1e6:.6f}"
                        if st["line_name"] == "other" else "")
                while True:
                    m = ask("Rest frequency (MHz)", prev)
                    if m == BACK:
                        break                  # back to the menu
                    try:
                        f = float(m) * 1e6
                    except ValueError:
                        print("  Not a number, try again.")
                        continue
                    if f <= 0:
                        print("  Must be positive, try again.")
                        continue
                    st["line_choice"] = v
                    st["line_name"] = "other"
                    st["rest_hz"] = f
                    return None
                continue
            print(f"  '{v}' is not valid: type a number from 1 to {n_other}.")

    def auto_names(st):
        date = st["selected"][0]["date_obs"][:10].replace("-", "")
        tag = st["target"].replace(" ", "_") if st["target"] else "session"
        return f"SALSA_{date}_{tag}", f"SALSA_{date}_{tag}.fro"

    def auto_dataset(st):
        glon = [i["glon"] for i in st["selected"]]
        glat = [i["glat"] for i in st["selected"]]
        return (f"GLon {min(glon):.1f} to {max(glon):.1f}, "
                f"GLat {min(glat):.1f} to {max(glat):.1f}")

    def text_step(key, prompt, auto=None):
        """Free-text question; with auto, the default follows the current
        selection until the operator types something different."""
        def step(st):
            if auto is not None and not st["custom"].get(key):
                st[key] = auto(st)
            v = ask(prompt, st[key])
            if v == BACK:
                return BACK
            if auto is not None:
                st["custom"][key] = (v != auto(st))
            elif v != st[key]:
                st["custom"][key] = True
            st[key] = v
            return None
        return step

    s = dict(
        selection=step_selection,
        line=step_line,
        target=text_step("target", "Target (optional, e.g. OH127.8+0.0)"),
        facility=text_step("facility", "Facility"),
        provider=text_step("provider", "Provider"),
        instrument=text_step("instrument", "Instrument"),
        survey=text_step("survey", "Survey"),
        dataset=text_step("dataset", "Dataset", auto_dataset),
        observation_id=text_step("observation_id", "Observation ID",
                                 lambda st: auto_names(st)[0]),
        output_filename=text_step("output_filename", "Output .fro filename",
                                  lambda st: auto_names(st)[1]),
        auto_dataset=auto_dataset, auto_names=auto_names,
    )
    return s


def run_steps(steps, st, start=0):
    i = start
    while i < len(steps):
        if steps[i](st) == BACK:
            if i == 0:
                print("  Already at the first question.")
            else:
                i -= 1
            continue
        i += 1


def main():
    input_dir = choose_input_dir()
    fits_files = sorted(glob.glob(os.path.join(input_dir, "*.fits")))
    if not fits_files:
        print(f"No FITS files found in {input_dir}.")
        return

    print(f"Found {len(fits_files)} FITS files in {input_dir}:\n")
    file_infos = [read_fits_header_info(fp) for fp in fits_files]
    for idx, info in enumerate(file_infos, start=1):
        date_display = (info["date_obs"].replace("T", " ")
                        if info["date_obs"] else "unknown")
        print(f"{idx:3d}. {os.path.basename(info['filepath'])}   "
              f"{date_display}   GLON={info['glon']:.1f}  "
              f"GLAT={info['glat']:.1f}  VELOSYS={info['velosys'] / 1e3:+.2f} km/s"
              f"  Az/El {'in header' if info['has_azel'] else 'MISSING'}")

    st = dict(selection=f"1-{len(file_infos)}", selected=None, axis=None, line_choice=None,
              line_name=None, rest_hz=np.nan, target="", facility="",
              provider=DEFAULT_PROVIDER, instrument="", survey="",
              dataset="", observation_id="", output_filename="", custom={})
    S = make_steps(file_infos)
    questions = [S["selection"], S["line"], S["target"], S["survey"],
                 S["output_filename"]]

    print("\nType < and Enter at any question to go back one step.\n")
    run_steps(questions, st)
    # Values not asked: fill them from the selection.
    for key in ("dataset", "observation_id", "output_filename"):
        if not st["custom"].get(key):
            st[key] = (S["auto_dataset"](st) if key == "dataset" else
                       S["auto_names"](st)[0 if key == "observation_id" else 1])

    review = [
        ("Selected files", lambda: f"{st['selection']} ({len(st['selected'])} files)",
         S["selection"]),
        ("Observed line", lambda: f"{st['line_name']}, rest frequency "
                                  f"{st['rest_hz'] / 1e6:.6f} MHz", S["line"]),
        ("Target", lambda: st["target"], S["target"]),
        ("Facility", lambda: st["facility"], S["facility"]),
        ("Provider", lambda: st["provider"], S["provider"]),
        ("Instrument", lambda: st["instrument"], S["instrument"]),
        ("Survey", lambda: st["survey"], S["survey"]),
        ("Dataset", lambda: st["dataset"], S["dataset"]),
        ("Observation ID", lambda: st["observation_id"], S["observation_id"]),
        ("Output filename", lambda: st["output_filename"], S["output_filename"]),
    ]

    while True:
        # Write the .fro next to the FITS files (the dataset folder) unless
        # a path was typed: that is where the HALO Launcher looks for it.
        out_path = st["output_filename"]
        if not os.path.dirname(out_path):
            out_path = os.path.join(input_dir, out_path)
        if os.path.exists(out_path):
            if not ask_yes_no(f"File '{out_path}' already exists. Overwrite?"):
                print("  File not overwritten: choose another filename.")
                S["output_filename"](st)
                continue

        print("\nPlease review the values before writing the file:")
        for k, (label, value, _) in enumerate(review, 1):
            print(f"  {k:2d}. {label + ':':17s}{value()}")
        print(f"      {'Written to:':17s}{out_path}")
        if ask_yes_no("Is everything correct?"):
            break
        while True:
            v = input(f"Number of the value to correct (1-{len(review)}, "
                      "Enter = show the review again): ").strip()
            if not v:
                break
            if v.isdigit() and 1 <= int(v) <= len(review):
                review[int(v) - 1][2](st)      # '<' here simply returns
                if int(v) == 1:                # new selection: refresh autos
                    for key in ("dataset", "observation_id", "output_filename"):
                        if not st["custom"].get(key):
                            st[key] = (S["auto_dataset"](st) if key == "dataset"
                                       else S["auto_names"](st)[
                                           0 if key == "observation_id" else 1])
                break
            print(f"  '{v}' is not valid: type a number from 1 to {len(review)}.")

    sel = st["selected"]
    ref = sel[0]
    axis = st["axis"]
    n = axis.size
    src_rest = ref["restfrq"]

    n_azel = sum(i["has_azel"] for i in sel)
    if n_azel == len(sel):
        pointing_source = ("commanded (SALSA FITS AZIMUTH/ELEVATIO, start of "
                           "integration)")
        azel_note = ("Az/El: commanded values from the FITS keywords "
                     "AZIMUTH/ELEVATIO (start of integration).")
    elif n_azel == 0:
        pointing_source = ("computed (Astropy, GLON/GLAT + DATE-OBS at Onsala; "
                           "AZIMUTH/ELEVATIO not in these FITS)")
        azel_note = ("Az/El: computed with Astropy from GLON/GLAT + DATE-OBS "
                     "at Onsala (the FITS have no AZIMUTH/ELEVATIO), not "
                     "telescope values.")
    else:
        pointing_source = (f"mixed: commanded (FITS AZIMUTH/ELEVATIO) for "
                           f"{n_azel} spectra, computed (Astropy) for "
                           f"{len(sel) - n_azel}")
        azel_note = (f"Az/El: commanded values from the FITS (AZIMUTH/ELEVATIO) "
                     f"for {n_azel} spectra; computed with Astropy from "
                     f"GLON/GLAT + DATE-OBS for the {len(sel) - n_azel} files "
                     "without them.")
    specs = sorted({i["specsys"] or "not given" for i in sel})
    doppler_note = (
        "vlsr_correction_mps: SALSA FITS keyword VELOSYS (m/s), one value per "
        f"spectrum, copied unchanged; reference frame (SPECSYS): "
        f"{', '.join(specs)}. Convention as used for SALSA data: "
        "v_LSR = v_topocentric (radio definition) + VELOSYS. NaN if the "
        "keyword is missing.")
    notes = (
        f"Converted from {len(sel)} SALSA FITS files via salsa_fits_to_fro.py. "
        f"Observed line chosen at import: {st['line_name']}, "
        f"{st['rest_hz'] / 1e6:.6f} MHz (the FITS RESTFRQ, "
        f"{src_rest / 1e6:.6f} MHz, is kept in Observation/source_restfrq_hz). "
        f"{azel_note} VELOSYS of each spectrum in Doppler/vlsr_correction_mps; "
        "the same correction computed independently with Astropy in "
        "Doppler/vlsr_correction_computed_mps."
    )

    create_fro_structure(
        out_path,
        frequency_axis_hz=axis,
        facility=st["facility"],
        instrument=st["instrument"],
        observation_id=st["observation_id"],
        survey=st["survey"],
        dataset=st["dataset"],
        provider=st["provider"],
        import_date=datetime.now(timezone.utc).isoformat(),
        observatory_name=f"Onsala Space Observatory - {st['instrument']}",
        latitude_deg=ONSALA_LOCATION.lat.deg,
        longitude_deg=ONSALA_LOCATION.lon.deg,
        elevation_m=ONSALA_LOCATION.height.to(u.m).value,
        receiver="SALSA remote receiver",
        antenna=f"{st['instrument']} (2.3m dish)",
        lna="n/a (remote facility, chain not controlled by HALO)",
        center_frequency_hz=float(axis[n // 2]),   # band centre channel
        bandwidth_hz=abs(n * float(ref["cdelt1"])),
        integration_time_s=np.nan,  # varies per spectrum, see Time/integration_time_s
        notes=notes,
        processed_by="HALO (salsa_fits_to_fro.py, FRO format v1.8)",
        pointing_source=pointing_source,
        target=st["target"],
        rest_frequency_hz=st["rest_hz"],
        source_restfrq_hz=src_rest,
        doppler_note=doppler_note,
        doppler_computed_note=DOPPLER_COMPUTED_NOTE,
        doppler_source_keyword="VELOSYS",
        overwrite=True,             # existence already confirmed above
    )

    print(f"\nCreated {out_path}, writing {len(sel)} spectra...\n")

    for info in sel:
        with fits.open(info["filepath"]) as hdul:
            data = hdul[0].data
        spectrum = data.reshape(-1).astype(np.float32)  # native byte order

        if info["has_azel"]:
            az_deg, el_deg, how = info["az"], info["el"], "commanded"
        else:
            az_deg, el_deg = compute_az_el(info["glon"], info["glat"],
                                           info["date_obs"])
            how = "computed"
        vcomp = compute_vlsr_correction(info["glon"], info["glat"],
                                        info["date_obs"])
        info["vcomp"] = vcomp
        append_raw_spectrum(
            out_path,
            spectrum,
            utc_time=datetime.fromisoformat(info["date_obs"]),
            az_deg=az_deg,
            elev_deg=el_deg,
            glon_deg=info["glon"],
            glat_deg=info["glat"],
            integration_time_s=info["obstime"],
            vlsr_correction_mps=info["velosys"],
            file_name=os.path.basename(info["filepath"]),
            vlsr_correction_computed_mps=vcomp,
        )
        print(f"  {os.path.basename(info['filepath'])}: "
              f"GLON={info['glon']:.1f} GLAT={info['glat']:.1f} -> "
              f"Az={az_deg:.2f} El={el_deg:.2f} ({how}), "
              f"integration {info['obstime']:.1f}s, "
              f"VELOSYS {info['velosys'] / 1e3:+.2f} km/s "
              f"(computed {vcomp / 1e3:+.2f})")

    bad = [i for i in sel if np.isfinite(i["velosys"])
           and abs(i["velosys"] - i["vcomp"]) > 1000.0]
    if bad:
        print(f"\nWARNING: VELOSYS differs from the computed LSR correction by "
              f"more than 1 km/s for {len(bad)} of {len(sel)} spectra "
              "(both values are stored):")
        for i in bad:
            print(f"  {os.path.basename(i['filepath'])}: GLON={i['glon']:.1f} "
                  f"GLAT={i['glat']:.1f}  VELOSYS {i['velosys'] / 1e3:+.2f}, "
                  f"computed {i['vcomp'] / 1e3:+.2f} km/s")

    print("\nDone.\n\n" + "=" * 42 + "\nSummary:\n" + "=" * 42)
    read_fro_summary(out_path)
    print("=" * 42)


if __name__ == "__main__":
    main()
