import h5py
import os
import numpy as np
from datetime import datetime, timezone
import textwrap


FORMAT_NAME = "FRO"
FORMAT_VERSION = "1.7"

# Sentinel for products with no defined observation epoch (mosaicked
# survey cubes such as HI4PI, THINGS, FEASTS, where each pixel combines
# multiple original observations). Pass utc_time=UNKNOWN_EPOCH explicitly;
# utc_time=None still means "not supplied" and defaults to the current time.
UNKNOWN_EPOCH = "unknown"

# FRO format v17 -- changes with respect to v16:
#   * Pointing: new per-row datasets ra_deg / dec_deg (measured equatorial
#     coordinates, available e.g. in SRT fitszilla files).
#   * New group "Quality": per-row flag_track (telescope tracking status,
#     0 = slewing / not on source) and flag_cal (calibration mark status).
#   * New group "Environment": per-row weather data (humidity, temperature,
#     pressure).
#   * New group "Calibration": per-RF-input parallel datasets
#     (polarization, cal_mark_temp_k, attenuation_db), fixed length equal
#     to the number of RF inputs, written once at creation time.
#   * Spectra: optional datasets cross_re / cross_im for the complex
#     cross-polarization product of full-Stokes backends (e.g. SARDARA).
#     NOTE: one row per integration, NOT per polarization, so they hold
#     half the rows of raw_spectrum when two polarizations are present.
#   * Observation: new attributes target, vlsr_mps, scan_id, subscan_id,
#     subscan_type, signal, schedule_name (scan metadata, useful for
#     position-switching pipelines).
# Design rule inherited from v16: row i of every per-row dataset
# (Spectra/raw_spectrum, Spectra/polarization, Time/*, Pointing/*,
# Quality/*, Environment/*) always refers to the same spectrum. With two
# polarizations each integration therefore produces two rows, and the
# shared quantities (time, pointing, flags, weather) are duplicated on
# both rows on purpose: redundancy is preferred over cross-indexing.


def now_utc():
    return datetime.now(timezone.utc)


def create_fro_v17_structure(
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
    # ----- new in v17 -----
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
        # New in v17: target and scan metadata.
        observation.attrs["target"] = target
        observation.attrs["vlsr_mps"] = vlsr_mps
        observation.attrs["scan_id"] = scan_id
        observation.attrs["subscan_id"] = subscan_id
        observation.attrs["subscan_type"] = subscan_type
        observation.attrs["signal"] = signal
        observation.attrs["schedule_name"] = schedule_name

        spectra = f.create_group("Spectra")
        spectra.create_dataset(
            "frequency_axis_hz",
            data=frequency_axis_hz,
            dtype="float64"
        )
        spectra.create_dataset(
            "raw_spectrum",
            shape=(0, channels),
            maxshape=(None, channels),
            chunks=(1, channels),
            dtype="float32"
        )
        spectra.create_dataset(
            "polarization",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype=string_dtype
        )
        # New in v17: optional complex cross-polarization product for
        # full-Stokes backends. One row per integration (see note above).
        if include_cross:
            spectra.create_dataset(
                "cross_re",
                shape=(0, channels),
                maxshape=(None, channels),
                chunks=(1, channels),
                dtype="float32"
            )
            spectra.create_dataset(
                "cross_im",
                shape=(0, channels),
                maxshape=(None, channels),
                chunks=(1, channels),
                dtype="float32"
            )
            spectra.attrs["cross_note"] = (
                "cross_re/cross_im hold one row per integration, "
                "not per polarization row of raw_spectrum."
            )

        time = f.create_group("Time")
        time.attrs["start_utc"] = ""
        time.attrs["end_utc"] = ""
        time.attrs["integration_time_s"] = integration_time_s
        time.create_dataset(
            "utc",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype=string_dtype
        )
        time.create_dataset(
            "unix_time",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        time.create_dataset(
            "integration_time_s",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )

        pointing = f.create_group("Pointing")
        pointing.attrs["pointing_source"] = pointing_source
        pointing.create_dataset(
            "azimuth_deg",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        pointing.create_dataset(
            "elevation_deg",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        pointing.create_dataset(
            "glon_deg",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        pointing.create_dataset(
            "glat_deg",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        # New in v17: measured equatorial coordinates (J2000), per row.
        pointing.create_dataset(
            "ra_deg",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        pointing.create_dataset(
            "dec_deg",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )

        # New in v17: per-row data quality flags from the telescope.
        quality = f.create_group("Quality")
        quality.attrs["flag_track_note"] = (
            "flag_track: 1 = telescope tracking / on source, "
            "0 = slewing or otherwise off source. Rows with 0 are kept "
            "on purpose (no data discarded at import time); filtering "
            "is left to the analysis stage."
        )
        quality.create_dataset(
            "flag_track",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="int32"
        )
        quality.create_dataset(
            "flag_cal",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )

        # New in v17: per-row weather data.
        environment = f.create_group("Environment")
        environment.create_dataset(
            "humidity_percent",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        environment.create_dataset(
            "temperature_c",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )
        environment.create_dataset(
            "pressure_mbar",
            shape=(0,),
            maxshape=(None,),
            chunks=(1,),
            dtype="float64"
        )

        # New in v17: per-RF-input calibration data, fixed length,
        # written once at creation time.
        calibration = f.create_group("Calibration")
        if all(cal_given):
            calibration.create_dataset(
                "polarization",
                data=np.asarray(cal_polarizations, dtype=object),
                dtype=string_dtype
            )
            calibration.create_dataset(
                "cal_mark_temp_k",
                data=np.asarray(cal_mark_temp_k, dtype=np.float64),
                dtype="float64"
            )
            calibration.create_dataset(
                "attenuation_db",
                data=np.asarray(attenuation_db, dtype=np.float64),
                dtype="float64"
            )


def append_raw_spectrum(filename, raw_spectrum, utc_time=None, az_deg=None, elev_deg=None, glon_deg=None, glat_deg=None, integration_time_s=None, polarization=None, ra_deg=None, dec_deg=None, flag_track=None, flag_cal=None, humidity_percent=None, temperature_c=None, pressure_mbar=None):
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

    with h5py.File(filename, "a") as f:
        dset = f["Spectra/raw_spectrum"]
        pol_dset = f["Spectra/polarization"]
        utc_dset = f["Time/utc"]
        unix_dset = f["Time/unix_time"]
        az_dset = f["Pointing/azimuth_deg"]
        elev_dset = f["Pointing/elevation_deg"]
        glon_dset = f["Pointing/glon_deg"]
        glat_dset = f["Pointing/glat_deg"]
        ra_dset = f["Pointing/ra_deg"]
        dec_dset = f["Pointing/dec_deg"]
        intt_dset = f["Time/integration_time_s"]
        track_dset = f["Quality/flag_track"]
        fcal_dset = f["Quality/flag_cal"]
        hum_dset = f["Environment/humidity_percent"]
        temp_dset = f["Environment/temperature_c"]
        press_dset = f["Environment/pressure_mbar"]

        channels = dset.shape[1]

        if raw_spectrum.shape != (channels,):
            raise ValueError(
                f"raw_spectrum must have {channels} channels, "
                f"but has {raw_spectrum.size}"
            )

        n = dset.shape[0]

        dset.resize((n + 1, channels))
        pol_dset.resize((n + 1,))
        utc_dset.resize((n + 1,))
        unix_dset.resize((n + 1,))
        az_dset.resize((n + 1,))
        elev_dset.resize((n + 1,))
        glon_dset.resize((n + 1,))
        glat_dset.resize((n + 1,))
        ra_dset.resize((n + 1,))
        dec_dset.resize((n + 1,))
        intt_dset.resize((n + 1,))
        track_dset.resize((n + 1,))
        fcal_dset.resize((n + 1,))
        hum_dset.resize((n + 1,))
        temp_dset.resize((n + 1,))
        press_dset.resize((n + 1,))

        dset[n, :] = raw_spectrum
        pol_dset[n] = polarization if polarization is not None else "N/A"
        utc_dset[n] = utc_string
        unix_dset[n] = unix_time
        az_dset[n] = az_deg if az_deg is not None else np.nan
        elev_dset[n] = elev_deg if elev_deg is not None else np.nan
        glon_dset[n] = glon_deg if glon_deg is not None else np.nan
        glat_dset[n] = glat_deg if glat_deg is not None else np.nan
        ra_dset[n] = ra_deg if ra_deg is not None else np.nan
        dec_dset[n] = dec_deg if dec_deg is not None else np.nan
        intt_dset[n] = integration_time_s if integration_time_s is not None else np.nan
        # flag_track defaults to 1 (tracking) when the source format does
        # not provide it, so that legacy-style data is not flagged out.
        track_dset[n] = flag_track if flag_track is not None else 1
        fcal_dset[n] = flag_cal if flag_cal is not None else np.nan
        hum_dset[n] = humidity_percent if humidity_percent is not None else np.nan
        temp_dset[n] = temperature_c if temperature_c is not None else np.nan
        press_dset[n] = pressure_mbar if pressure_mbar is not None else np.nan

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
# the original standard ISO 8601 format; only this print-time
# representation is adjusted.
def format_utc_display(utc_string):
    if not utc_string:
        return utc_string
    base, _, offset = utc_string.partition("+")
    base = base.split(".")[0]
    return base.replace("T", " ") + ("+" + offset if offset else "")


def read_fro_summary(filename):
    label_width = 20
    with h5py.File(filename, "r") as f:
        raw = f["Spectra/raw_spectrum"]
        print(f"{'Format:':<{label_width}}{f.attrs['format_name']}")
        print(f"{'Version:':<{label_width}}{f.attrs['format_version']}")
        if "Source" in f:
            print(f"{'Facility:':<{label_width}}{f['Source'].attrs.get('facility', '')}")
            print(f"{'Instrument:':<{label_width}}{f['Source'].attrs.get('instrument', '')}")
            print(f"{'Survey:':<{label_width}}{f['Source'].attrs.get('survey', '')}")
            print(f"{'Dataset:':<{label_width}}{f['Source'].attrs.get('dataset', '')}")
            print(f"{'Observation ID:':<{label_width}}{f['Source'].attrs.get('observation_id', '')}")
            print(f"{'Provider:':<{label_width}}{f['Source'].attrs.get('provider', '')}")
            print(f"{'Import date:':<{label_width}}{format_utc_display(f['Source'].attrs.get('import_date', ''))}")
            print(f"{'Processed by:':<{label_width}}{f['Source'].attrs.get('processed_by', '')}")
        else:
            print("Source: not present")
        print(f"{'Observatory:':<{label_width}}{f['Observatory'].attrs['name']}")
        print(f"{'Receiver:':<{label_width}}{f['Hardware'].attrs['receiver']}")
        print(f"{'Antenna:':<{label_width}}{f['Hardware'].attrs['antenna']}")
        print(f"{'LNA:':<{label_width}}{f['Hardware'].attrs['lna']}")
        target = f["Observation"].attrs.get("target", "")
        if target:
            print(f"{'Target:':<{label_width}}{target}")
        scan_id = f["Observation"].attrs.get("scan_id", "")
        subscan_id = f["Observation"].attrs.get("subscan_id", "")
        subscan_type = f["Observation"].attrs.get("subscan_type", "")
        if scan_id or subscan_id or subscan_type:
            print(f"{'Scan:':<{label_width}}id={scan_id} subscan={subscan_id} type={subscan_type}")
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
        print(f"{'Temp. calibration:':<{label_width}}{f['Observation'].attrs.get('temperature_calibration_note', '')}")
        # Notes is a free-text field (e.g. explaining how pointing was
        # derived, or listing source files) and can be long. Word-wrap it
        # to a fixed width and indent continuation lines under the same
        # label column used by every other field above, for readability.
        if "Header" in f and f["Header"].attrs.get("notes", ""):
            wrapped = textwrap.wrap(f['Header'].attrs['notes'], width=80)
            print(f"{'Notes:':<{label_width}}{wrapped[0]}")
            for line in wrapped[1:]:
                print(f"{'':<{label_width}}{line}")