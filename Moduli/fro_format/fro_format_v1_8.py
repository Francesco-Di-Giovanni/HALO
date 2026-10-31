"""
fro_format_v1_8.py - FRO format, version 1.8 (HDF5).

Superset of version 1.7: every v1.7 field and argument is unchanged; the
additions are optional and default to NaN / empty when not supplied.

Changes with respect to v1.7:
  * Observation: new attributes
      rest_frequency_hz  - rest frequency of the line actually observed,
                           chosen at import time (e.g. OH 1612.231 MHz);
      source_restfrq_hz  - rest frequency written in the source files,
                           kept as found even when it does not match the
                           observed line (e.g. SALSA always writes HI).
  * New group "Doppler": per-row dataset vlsr_correction_mps, the
    line-of-sight correction from the topocentric to the LSR frame as
    supplied by the source data (NaN if not supplied). The group attribute
    "vlsr_correction_note" states where the value comes from and its sign
    convention. Per-row vlsr_correction_computed_mps holds the same
    correction computed independently (method in the attribute
    "vlsr_computed_note"; NaN if not computed), so that wrong source
    values can be detected; the source value is never replaced.
  * Source: new per-row dataset file_name (name of the source file of each
    row), for traceability.

Function names carry no version number from v1.8 on (create_fro_structure),
so that a future version only requires changing the import line.

Design rule inherited from v16/v17: row i of every per-row dataset
(Spectra/raw_spectrum, Spectra/polarization, Time/*, Pointing/*, Quality/*,
Environment/*, Doppler/*, Source/file_name) always refers to the same
spectrum. With two polarizations each integration produces two rows and the
shared quantities are duplicated on both rows on purpose: redundancy is
preferred over cross-indexing.

v1.7 notes, still valid:
  * Pointing: per-row ra_deg / dec_deg (measured equatorial coordinates).
  * Group "Quality": per-row flag_track (0 = slewing / not on source) and
    flag_cal (calibration mark status).
  * Group "Environment": per-row humidity, temperature, pressure.
  * Group "Calibration": per-RF-input datasets (polarization,
    cal_mark_temp_k, attenuation_db), written once at creation time.
  * Spectra: optional cross_re / cross_im (one row per integration, NOT per
    polarization row of raw_spectrum).
  * Observation: target, vlsr_mps (source velocity, one value per file),
    scan_id, subscan_id, subscan_type, signal, schedule_name.
"""
import h5py
import os
import numpy as np
from datetime import datetime, timezone
import textwrap


FORMAT_NAME = "FRO"
FORMAT_VERSION = "1.8"

# Sentinel for products with no defined observation epoch (mosaicked
# survey cubes such as HI4PI, THINGS, FEASTS, where each pixel combines
# multiple original observations). Pass utc_time=UNKNOWN_EPOCH explicitly;
# utc_time=None still means "not supplied" and defaults to the current time.
UNKNOWN_EPOCH = "unknown"

DEFAULT_DOPPLER_NOTE = (
    "vlsr_correction_mps: line-of-sight velocity correction from the "
    "topocentric to the LSR frame, per row, as supplied by the source data "
    "(NaN if not supplied). Convention: v_LSR = v_topocentric (radio "
    "definition) + vlsr_correction."
)


def now_utc():
    return datetime.now(timezone.utc)


def _row_dataset(group, name, dtype):
    """Create an empty, extendable per-row dataset."""
    group.create_dataset(name, shape=(0,), maxshape=(None,), chunks=(1,),
                         dtype=dtype)


def create_fro_structure(
    filename,
    frequency_axis_hz,
    facility,
    instrument,
    observation_id,
    survey="",
    dataset="",
    provider="",
    import_date="",
    observatory_name="Francesco Radio Observatory",
    latitude_deg=np.nan,
    longitude_deg=np.nan,
    elevation_m=np.nan,
    receiver="unknown",
    antenna="unknown",
    lna="unknown",
    center_frequency_hz=np.nan,
    sample_rate_hz=np.nan,
    bandwidth_hz=np.nan,
    integration_time_s=np.nan,
    notes="",
    temperature_calibration_note="Antenna temperature (K) based on assumed Tsys ~300K; not an absolute calibration.",
    processed_by="",
    pointing_source="unknown",
    target="",
    vlsr_mps=np.nan,
    scan_id="",
    subscan_id="",
    subscan_type="",
    signal="",
    schedule_name="",
    cal_polarizations=None,
    cal_mark_temp_k=None,
    attenuation_db=None,
    include_cross=False,
    # ----- new in v1.8 -----
    rest_frequency_hz=np.nan,
    source_restfrq_hz=np.nan,
    doppler_note=DEFAULT_DOPPLER_NOTE,
    doppler_computed_note="",
    doppler_source_keyword="",
    overwrite=False
):
    frequency_axis_hz = np.asarray(frequency_axis_hz, dtype=np.float64)
    channels = frequency_axis_hz.size

    if channels < 1:
        raise ValueError("frequency_axis_hz must not be empty")
    if not facility:
        raise ValueError("facility is mandatory")
    if not instrument:
        raise ValueError("instrument is mandatory")
    if not observation_id:
        raise ValueError("observation_id is mandatory")

    # Calibration inputs must be either all absent or all present with
    # matching lengths (one entry per RF input, e.g. HLP and VLP).
    cal_given = [x is not None for x in
                 (cal_polarizations, cal_mark_temp_k, attenuation_db)]
    if any(cal_given) and not all(cal_given):
        raise ValueError(
            "cal_polarizations, cal_mark_temp_k and attenuation_db "
            "must be provided together"
        )
    if all(cal_given):
        n_inputs = len(cal_polarizations)
        if not (len(cal_mark_temp_k) == n_inputs
                and len(attenuation_db) == n_inputs):
            raise ValueError(
                "cal_polarizations, cal_mark_temp_k and attenuation_db "
                "must have the same length"
            )

    created = now_utc().isoformat()
    string_dtype = h5py.string_dtype(encoding="utf-8")
    if os.path.exists(filename) and not overwrite:
        raise FileExistsError(f"Output file already exists: {filename}")

    with h5py.File(filename, "w") as f:
        f.attrs["format_name"] = FORMAT_NAME
        f.attrs["format_version"] = FORMAT_VERSION
        f.attrs["created_utc"] = created

        header = f.create_group("Header")
        header.attrs["format_name"] = FORMAT_NAME
        header.attrs["format_version"] = FORMAT_VERSION
        header.attrs["created_utc"] = created
        header.attrs["notes"] = notes

        source = f.create_group("Source")
        source.attrs["facility"] = facility
        source.attrs["instrument"] = instrument
        source.attrs["survey"] = survey
        source.attrs["dataset"] = dataset
        source.attrs["observation_id"] = observation_id
        source.attrs["provider"] = provider
        source.attrs["import_date"] = import_date
        source.attrs["processed_by"] = processed_by
        # New in v1.8: source file of each row.
        _row_dataset(source, "file_name", string_dtype)

        observatory = f.create_group("Observatory")
        observatory.attrs["name"] = observatory_name
        observatory.attrs["latitude_deg"] = latitude_deg
        observatory.attrs["longitude_deg"] = longitude_deg
        observatory.attrs["elevation_m"] = elevation_m

        hardware = f.create_group("Hardware")
        hardware.attrs["receiver"] = receiver
        hardware.attrs["antenna"] = antenna
        hardware.attrs["lna"] = lna

        observation = f.create_group("Observation")
        observation.attrs["center_frequency_hz"] = center_frequency_hz
        observation.attrs["sample_rate_hz"] = sample_rate_hz
        observation.attrs["bandwidth_hz"] = bandwidth_hz
        observation.attrs["channels"] = channels
        observation.attrs["integration_time_s"] = integration_time_s
        observation.attrs["temperature_calibration_note"] = temperature_calibration_note
        observation.attrs["target"] = target
        observation.attrs["vlsr_mps"] = vlsr_mps
        observation.attrs["scan_id"] = scan_id
        observation.attrs["subscan_id"] = subscan_id
        observation.attrs["subscan_type"] = subscan_type
        observation.attrs["signal"] = signal
        observation.attrs["schedule_name"] = schedule_name
        # New in v1.8: observed line and the value found in the source.
        observation.attrs["rest_frequency_hz"] = rest_frequency_hz
        observation.attrs["source_restfrq_hz"] = source_restfrq_hz

        spectra = f.create_group("Spectra")
        spectra.create_dataset("frequency_axis_hz", data=frequency_axis_hz,
                               dtype="float64")
        spectra.create_dataset("raw_spectrum", shape=(0, channels),
                               maxshape=(None, channels),
                               chunks=(1, channels), dtype="float32")
        _row_dataset(spectra, "polarization", string_dtype)
        if include_cross:
            for name in ("cross_re", "cross_im"):
                spectra.create_dataset(name, shape=(0, channels),
                                       maxshape=(None, channels),
                                       chunks=(1, channels), dtype="float32")
            spectra.attrs["cross_note"] = (
                "cross_re/cross_im hold one row per integration, "
                "not per polarization row of raw_spectrum."
            )

        time = f.create_group("Time")
        time.attrs["start_utc"] = ""
        time.attrs["end_utc"] = ""
        time.attrs["integration_time_s"] = integration_time_s
        _row_dataset(time, "utc", string_dtype)
        _row_dataset(time, "unix_time", "float64")
        _row_dataset(time, "integration_time_s", "float64")

        pointing = f.create_group("Pointing")
        pointing.attrs["pointing_source"] = pointing_source
        for name in ("azimuth_deg", "elevation_deg", "glon_deg", "glat_deg",
                     "ra_deg", "dec_deg"):
            _row_dataset(pointing, name, "float64")

        quality = f.create_group("Quality")
        quality.attrs["flag_track_note"] = (
            "flag_track: 1 = telescope tracking / on source, "
            "0 = slewing or otherwise off source. Rows with 0 are kept "
            "on purpose (no data discarded at import time); filtering "
            "is left to the analysis stage."
        )
        _row_dataset(quality, "flag_track", "int32")
        _row_dataset(quality, "flag_cal", "float64")

        environment = f.create_group("Environment")
        for name in ("humidity_percent", "temperature_c", "pressure_mbar"):
            _row_dataset(environment, name, "float64")

        # New in v1.8: per-row LSR correction.
        doppler = f.create_group("Doppler")
        doppler.attrs["vlsr_correction_note"] = doppler_note
        _row_dataset(doppler, "vlsr_correction_mps", "float64")
        doppler.attrs["vlsr_computed_note"] = doppler_computed_note
        doppler.attrs["source_keyword"] = doppler_source_keyword
        _row_dataset(doppler, "vlsr_correction_computed_mps", "float64")

        calibration = f.create_group("Calibration")
        if all(cal_given):
            calibration.create_dataset(
                "polarization",
                data=np.asarray(cal_polarizations, dtype=object),
                dtype=string_dtype)
            calibration.create_dataset(
                "cal_mark_temp_k",
                data=np.asarray(cal_mark_temp_k, dtype=np.float64),
                dtype="float64")
            calibration.create_dataset(
                "attenuation_db",
                data=np.asarray(attenuation_db, dtype=np.float64),
                dtype="float64")


def append_raw_spectrum(filename, raw_spectrum, utc_time=None, az_deg=None,
                        elev_deg=None, glon_deg=None, glat_deg=None,
                        integration_time_s=None, polarization=None,
                        ra_deg=None, dec_deg=None, flag_track=None,
                        flag_cal=None, humidity_percent=None,
                        temperature_c=None, pressure_mbar=None,
                        vlsr_correction_mps=None, file_name=None,
                        vlsr_correction_computed_mps=None):
    raw_spectrum = np.asarray(raw_spectrum, dtype=np.float32)

    if utc_time is None:
        utc_time = now_utc()
    if isinstance(utc_time, str):
        if utc_time != UNKNOWN_EPOCH:
            raise ValueError(
                f"utc_time as a string is only allowed as the explicit "
                f"sentinel {UNKNOWN_EPOCH!r}; pass a datetime object otherwise"
            )
        utc_string = UNKNOWN_EPOCH
        unix_time = np.nan
    else:
        utc_string = utc_time.isoformat()
        unix_time = utc_time.timestamp()

    def nan_if_none(x):
        return x if x is not None else np.nan

    # Per-row values, keyed by the dataset they go to.
    row = {
        "Spectra/polarization": polarization if polarization is not None else "N/A",
        "Time/utc": utc_string,
        "Time/unix_time": unix_time,
        "Time/integration_time_s": nan_if_none(integration_time_s),
        "Pointing/azimuth_deg": nan_if_none(az_deg),
        "Pointing/elevation_deg": nan_if_none(elev_deg),
        "Pointing/glon_deg": nan_if_none(glon_deg),
        "Pointing/glat_deg": nan_if_none(glat_deg),
        "Pointing/ra_deg": nan_if_none(ra_deg),
        "Pointing/dec_deg": nan_if_none(dec_deg),
        # flag_track defaults to 1 (tracking) when the source format does
        # not provide it, so that legacy-style data is not flagged out.
        "Quality/flag_track": flag_track if flag_track is not None else 1,
        "Quality/flag_cal": nan_if_none(flag_cal),
        "Environment/humidity_percent": nan_if_none(humidity_percent),
        "Environment/temperature_c": nan_if_none(temperature_c),
        "Environment/pressure_mbar": nan_if_none(pressure_mbar),
        "Doppler/vlsr_correction_mps": nan_if_none(vlsr_correction_mps),
        "Doppler/vlsr_correction_computed_mps": nan_if_none(vlsr_correction_computed_mps),
        "Source/file_name": file_name if file_name is not None else "",
    }

    with h5py.File(filename, "a") as f:
        dset = f["Spectra/raw_spectrum"]
        channels = dset.shape[1]
        if raw_spectrum.shape != (channels,):
            raise ValueError(
                f"raw_spectrum must have {channels} channels, "
                f"but has {raw_spectrum.size}"
            )
        n = dset.shape[0]
        dset.resize((n + 1, channels))
        dset[n, :] = raw_spectrum
        for path, value in row.items():
            d = f[path]
            d.resize((n + 1,))
            d[n] = value

        if n == 0:
            f["Time"].attrs["start_utc"] = utc_string
        f["Time"].attrs["end_utc"] = utc_string


def append_cross_spectrum(filename, cross_re, cross_im):
    # Appends one integration of the complex cross-polarization product.
    # Only valid for files created with include_cross=True.
    cross_re = np.asarray(cross_re, dtype=np.float32)
    cross_im = np.asarray(cross_im, dtype=np.float32)

    with h5py.File(filename, "a") as f:
        if "Spectra/cross_re" not in f:
            raise ValueError(
                "this FRO file was created without cross datasets "
                "(include_cross=False)"
            )
        re_dset = f["Spectra/cross_re"]
        im_dset = f["Spectra/cross_im"]
        channels = re_dset.shape[1]
        if cross_re.shape != (channels,) or cross_im.shape != (channels,):
            raise ValueError(
                f"cross_re and cross_im must have {channels} channels"
            )
        n = re_dset.shape[0]
        re_dset.resize((n + 1, channels))
        im_dset.resize((n + 1, channels))
        re_dset[n, :] = cross_re
        im_dset[n, :] = cross_im


# Helper function: format an ISO 8601 UTC string for display only.
# Strips microseconds and replaces the "T" separator with a space,
# purely for readability on screen. The stored data always keeps
# the original standard ISO 8601 format.
def format_utc_display(utc_string):
    if not utc_string:
        return utc_string
    base, _, offset = utc_string.partition("+")
    base = base.split(".")[0]
    return base.replace("T", " ") + ("+" + offset if offset else "")


def read_fro_summary(filename):
    """Print a summary; works on v1.8 files and on older versions (fields
    that do not exist in the file are simply not shown)."""
    label_width = 27
    with h5py.File(filename, "r") as f:
        raw = f["Spectra/raw_spectrum"]
        print(f"{'Format:':<{label_width}}{f.attrs['format_name']}")
        print(f"{'Version:':<{label_width}}{f.attrs['format_version']}")
        if "Source" in f:
            src = f["Source"].attrs
            print(f"{'Facility:':<{label_width}}{src.get('facility', '')}")
            print(f"{'Instrument:':<{label_width}}{src.get('instrument', '')}")
            print(f"{'Survey:':<{label_width}}{src.get('survey', '')}")
            print(f"{'Dataset:':<{label_width}}{src.get('dataset', '')}")
            print(f"{'Observation ID:':<{label_width}}{src.get('observation_id', '')}")
            print(f"{'Provider:':<{label_width}}{src.get('provider', '')}")
            print(f"{'Import date:':<{label_width}}{format_utc_display(src.get('import_date', ''))}")
            print(f"{'Processed by:':<{label_width}}{src.get('processed_by', '')}")
        else:
            print("Source: not present")
        print(f"{'Observatory:':<{label_width}}{f['Observatory'].attrs['name']}")
        print(f"{'Receiver:':<{label_width}}{f['Hardware'].attrs['receiver']}")
        print(f"{'Antenna:':<{label_width}}{f['Hardware'].attrs['antenna']}")
        print(f"{'LNA:':<{label_width}}{f['Hardware'].attrs['lna']}")
        obs = f["Observation"].attrs
        target = obs.get("target", "")
        if target:
            print(f"{'Target:':<{label_width}}{target}")
        scan_id = obs.get("scan_id", "")
        subscan_id = obs.get("subscan_id", "")
        subscan_type = obs.get("subscan_type", "")
        if scan_id or subscan_id or subscan_type:
            print(f"{'Scan:':<{label_width}}id={scan_id} subscan={subscan_id} type={subscan_type}")
        rest = float(obs.get("rest_frequency_hz", np.nan))
        src_rest = float(obs.get("source_restfrq_hz", np.nan))
        if np.isfinite(rest):
            print(f"{'Rest frequency:':<{label_width}}{rest / 1e6:.6f} MHz")
        if np.isfinite(src_rest) and src_rest != rest:
            print(f"{'Source RESTFRQ:':<{label_width}}{src_rest / 1e6:.6f} MHz "
                  "(as found in the source files)")
        print(f"{'Spectra:':<{label_width}}{raw.shape[0]}")
        print(f"{'Channels:':<{label_width}}{raw.shape[1]}")
        if "Spectra" in f and "cross_re" in f["Spectra"]:
            print(f"{'Cross rows:':<{label_width}}{f['Spectra/cross_re'].shape[0]}")
        print(f"{'Start UTC:':<{label_width}}{format_utc_display(f['Time'].attrs['start_utc'])}")
        print(f"{'End UTC:':<{label_width}}{format_utc_display(f['Time'].attrs['end_utc'])}")
        if "Pointing" in f:
            print(f"{'Pointing source:':<{label_width}}{f['Pointing'].attrs.get('pointing_source', 'unknown')}")
        if raw.shape[0] > 0 and "Pointing" in f and "glon_deg" in f["Pointing"]:
            glon = f["Pointing/glon_deg"][:]
            glat = f["Pointing/glat_deg"][:]
            if not np.all(np.isnan(glon)):
                print(f"{'GLON range:':<{label_width}}{np.nanmin(glon)} to {np.nanmax(glon)}")
                print(f"{'GLAT range:':<{label_width}}{np.nanmin(glat)} to {np.nanmax(glat)}")
        if raw.shape[0] > 0 and "Pointing" in f and "ra_deg" in f["Pointing"]:
            ra = f["Pointing/ra_deg"][:]
            dec = f["Pointing/dec_deg"][:]
            if not np.all(np.isnan(ra)):
                print(f"{'RA range:':<{label_width}}{np.nanmin(ra):.4f} to {np.nanmax(ra):.4f}")
                print(f"{'DEC range:':<{label_width}}{np.nanmin(dec):.4f} to {np.nanmax(dec):.4f}")
        if raw.shape[0] > 0 and "Doppler" in f and "vlsr_correction_mps" in f["Doppler"]:
            dop = f["Doppler/vlsr_correction_mps"][:]
            kw = str(f["Doppler"].attrs.get("source_keyword", "")) or "source"
            lab_c = f"LSR correction ({kw}):"
            lab_k = f"LSR check ({kw}):"
            if np.all(np.isnan(dop)):
                print(f"{lab_c:<{label_width}}not supplied")
            else:
                print(f"{lab_c:<{label_width}}{np.nanmin(dop) / 1e3:+.2f} "
                      f"to {np.nanmax(dop) / 1e3:+.2f} km/s")
            if "vlsr_correction_computed_mps" in f["Doppler"]:
                cmp_ = f["Doppler/vlsr_correction_computed_mps"][:]
                if not np.all(np.isnan(cmp_)):
                    diff = (dop - cmp_) / 1e3
                    nbad = int(np.sum(np.abs(diff) > 1.0))
                    print(f"{lab_k:<{label_width}}{kw} - computed: median "
                          f"{np.nanmedian(diff):+.2f} km/s, max |diff| "
                          f"{np.nanmax(np.abs(diff)):.2f} km/s, {nbad} of {diff.size} "
                          "spectra above 1 km/s")
        if raw.shape[0] > 0 and "Quality" in f and "flag_track" in f["Quality"]:
            track = f["Quality/flag_track"][:]
            n_off = int(np.sum(track == 0))
            print(f"{'Rows not tracking:':<{label_width}}{n_off} of {track.size}")
        if "Calibration" in f and "polarization" in f["Calibration"]:
            pols = [p.decode() if isinstance(p, bytes) else p for p in f["Calibration/polarization"][:]]
            temps = f["Calibration/cal_mark_temp_k"][:]
            pairs = ", ".join(f"{p}={t:g}K" for p, t in zip(pols, temps))
            print(f"{'Cal marks:':<{label_width}}{pairs}")
        if raw.shape[0] > 0 and "Spectra" in f and "polarization" in f["Spectra"]:
            pols = [p.decode() if isinstance(p, bytes) else p for p in f["Spectra/polarization"][:]]
            unique_pols = sorted(set(pols))
            print(f"{'Polarizations:':<{label_width}}{', '.join(unique_pols)}")
        if raw.shape[0] > 0 and "Source" in f and "file_name" in f["Source"]:
            names = [x.decode() if isinstance(x, bytes) else x for x in f["Source/file_name"][:]]
            named = [x for x in names if x]
            if named:
                print(f"{'Source files:':<{label_width}}{len(set(named))} "
                      f"({named[0]} ... {named[-1]})")
        print(f"{'Temp. calibration:':<{label_width}}{obs.get('temperature_calibration_note', '')}")
        # Notes is a free-text field and can be long: word-wrap it under
        # the same label column used by every other field above.
        if "Header" in f and f["Header"].attrs.get("notes", ""):
            wrapped = textwrap.wrap(f['Header'].attrs['notes'], width=80)
            print(f"{'Notes:':<{label_width}}{wrapped[0]}")
            for line in wrapped[1:]:
                print(f"{'':<{label_width}}{line}")
