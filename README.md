# H.A.L.O. — Hydrogen Atomic Line Observatory

**Bolzano/Bozen, Italy — 46.4946°N, 11.3353°E, 334 m a.s.l.**

---

### What is HALO?

HALO is an amateur radio astronomy project based in Bolzano/Bozen, Italy, with an unusual ambition: to build a **telescope-agnostic, format-agnostic spectral data pipeline** that can ingest observations from any radio telescope — from a 1.2-metre 3D-printed dish in a backyard to the 500-metre FAST — the largest single-dish radio telescope ever built — and the 27-antenna VLA — and archive them in a single, unified HDF5 format for scientific analysis.

The project started in January 2026 with a simple 3D-printed parabolic antenna and an Airspy R2 SDR receiver, targeting the HI 21-cm line. It has since grown far beyond its origins.

---

### From a 3D-printed dish to Effelsberg

The journey so far:

- **January 2026** — Project begins. A 1.2-metre parabolic dish, 3D-printed and hand-assembled, connected to an Airspy R2 SDR and a Nooelec SawBird H1 LNA. First HI spectra acquired with the ezRA suite (Ted Cline, BAA/SARA).
- **Early 2026** — First contact with **SALSA** (Onsala Space Observatory, Sweden). Remote observations of the Galactic plane at 1.4 GHz. A bug in the SALSA slewing/integration system, discovered during these sessions, was reported and fixed by the Onsala team within hours (v1.1.8). The local SDR receiver configuration (8192 channels @ 2.5 MSps, Airspy R2) was validated during these sessions, achieving a spectral velocity resolution of **64 m/s per channel** at 1.4 GHz — a figure comparable to many professional survey instruments, obtained with consumer-grade hardware and long integration times.
- **Mid 2026** — The FRO format (HDF5 v17) is defined. Parsers written for SALSA, **GBT** (Green Bank Telescope, 100m), **FAST/FEASTS** (Five-hundred-metre Aperture Spherical Telescope — at 500 metres in diameter, the largest single-dish radio telescope ever built), **HI4PI**, **LAB**, **EBHIS**, **Parkes/GASS**, **ATCA**, **VLASS** (VLA Sky Survey), and **JCMT/COHRS** (12CO J=3-2 at 345.796 GHz — the first molecular line integrated into the pipeline, demonstrating that HALO is not limited to the HI line).
- **August 2026** — First contact with **Effelsberg** (Max Planck Institute for Radio Astronomy, Bonn). Real MBFITS HI data of **Holmberg 1** provided by Dr. Uwe Bach (MPIfR) for parser development and validation.

---

### Philosophy

HALO is built on a few core principles:

- **No hardcoded assumptions.** Every parameter — frequency axis, bandwidth, number of channels, coordinates, epoch — is read dynamically from the source file. Nothing is assumed, everything is verified.
- **Honest metadata.** If a quantity is not measured or not defined (e.g. the epoch of a mosaic product), it is explicitly flagged as unknown rather than filled with a plausible-looking value.
- **Raw data only.** The `.fro` archive format stores only the raw, uncalibrated spectra as delivered by the telescope. Smoothing, baseline subtraction, and calibration are reserved for the analysis stage.
- **No format lock-in.** The FRO HDF5 format is the hub, not the destination. Export converters to ezRA, SDFITS, and other formats are part of the roadmap.
- **Open to any spectral line.** Despite the "HALO" name, the pipeline handles any spectral line at any frequency — from HI at 1.4 GHz to 12CO(3-2) at 345.796 GHz and beyond.

---

### Telescopes and surveys currently supported

| Telescope / Survey | Frequency | Line | Parser |
|---|---|---|---|
| Local (Airspy R2 + 1.2m dish) | 1.4 GHz | HI | ezRA / ezColAirspy |
| SALSA (Onsala, 2×2.3m) | 1.4 GHz | HI | `salsa_fits_to_fro.py` |
| FAST/FEASTS (500m) | 1.4 GHz | HI | `feasts_cube_to_fro.py` |
| GBT (100m, Green Bank) | 1.4 GHz | HI, H₂O | `gbt_sdfits_to_fro.py` |
| Effelsberg (100m, Bonn) | 1.4 GHz | HI | `mbfits_to_fro.py` *(in dev)* |
| Parkes/GASS (64m) | 1.4 GHz | HI | `parkes_gass_to_fro.py` |
| HI4PI (all-sky survey) | 1.4 GHz | HI | `hi4pi_to_fro.py` |
| LAB Survey | 1.4 GHz | HI | `lab_to_fro.py` |
| EBHIS | 1.4 GHz | HI | `ebhis_to_fro.py` |
| ATCA (archival) | various | various | *(in dev)* |
| VLA/VLASS | 3 GHz | continuum | `view_vlass.py` |
| JCMT/COHRS (15m, Maunakea) | 345.796 GHz | ¹²CO(3-2) | `cohrs_to_fro.py` |

---

### Interactive viewers

| Viewer | Description | Link |
|---|---|---|
| FASHI DR2 — 3D HI Universe | 156,411 extragalactic HI sources from FAST, rendered in 3D. Colour-coded by HI mass. Rotate, zoom, and hover for source details. | [Open viewer](https://htmlpreview.github.io/?https://github.com/Francesco-Di-Giovanni/HALO/blob/main/FASHI_DR2_3D_universe.html) |

---

### The FRO HDF5 format (v17)

The core of the pipeline is `fro_format_v17.py`, which defines a common HDF5 structure for storing spectral observations from any telescope. Key groups:

- **Spectra** — raw spectra, frequency axis in Hz, optional cross-polarisation terms
- **Pointing** — Az/El, RA/Dec, GLon/GLat per integration
- **Time** — UTC timestamps, Unix time, integration duration
- **Source** — provenance metadata (facility, instrument, survey, dataset, provider)
- **Observatory** — physical location of the telescope
- **Quality** — tracking and calibration flags
- **Environment** — temperature, humidity, pressure per integration
- **Calibration** — per-feed polarisation and calibration parameters
- **Observation** — receiver, centre frequency, calibration notes
- **Notes** — free-text provenance and processing notes

The format uses the `UNKNOWN_EPOCH` sentinel for mosaic products where a per-pixel epoch is undefined (e.g. HI4PI, COHRS), ensuring metadata honesty at all times.

---

### Software and hardware

**Local receiver:**
- Airspy R2 SDR (2.5 / 10 MSps)
- Nooelec SawBird H1 LNA
- 1.2-metre parabolic dish (3D-printed)
- Fujitsu ESPRIMO Q958 mini-PC, Ubuntu, Python 3.14

**Analysis:**
- [ezRA suite](https://github.com/tedcline/ezRA) — primary analysis and visualisation tool, developed by Ted Cline (N0RQV)
- [DSPIRA](https://github.com/WVURAIL/gr-radio_astro) (Digital Signal Processing in Radio Astronomy) — GNU Radio flowgraph framework developed by West Virginia University (WVURAIL); used as the basis for the local SDR spectrometer, inspired by Andrew Sutkowski's GNU Radio interface for integration with ezRA. Thanks to Dr. Andrew Thornett (BAA/SARA) for his Monday online meetings.
- Custom viewers for each telescope (spectrum plots, l-v diagrams, moment maps, peak temperature maps)

**Roadmap — Italian radio telescopes:**
Parsers for the three major Italian radio telescopes are planned:
- **SRT** (Sardinia Radio Telescope, 64m, INAF)
- **Noto** (32m, INAF, Sicily)
- **Medicina** (32m, INAF, Bologna)

---

### Collaborators and acknowledgements

- **Ted Cline, N0RQV** (Little Thompson Observatory, Berthoud, CO-USA) — author of the ezRA suite; ongoing collaboration on parser development and GBT data
- **Dr. Eskil Varenius** (Onsala Space Observatory) — SALSA support and observing methodology
- **Dr. Uwe Bach** (MPIfR, Effelsberg) — provided real Effelsberg MBFITS HI data of Holmberg 1 for parser development
- **Dr. Chuan-Peng Zhang** (NAOC) — provided FAST M33 HI coordinate data for parser development and validation
- **Dr. Jing Wang** (KIAA, Peking University; FEASTS team lead) — provided permission to use FEASTS HI data of NGC 628 in the HALO pipeline and repository
- **Dr. Mario Sandri** (Unione Astrofili Italiani; Associazione Italiana di Fisica; Phoenix APS) — astrophysicist; one of the inspirations behind this project, first encountered at ICARA 2025 (Italian Congress of Amateur Radio Astronomy, Pordenone, October 2025)
- **Dr. Andrew Thornett** (BAA, SARA) Lichfield Radio Observatory (LRO), Lichfield, UK
- **Phoenix APS** (Cles, TN - Italy) — Amateur astronomy, astrophotography, radio astronomy and radio amateur club
- **Claude AI** (Anthropic) — developed in collaboration with Claude AI; all parsers, format specification, viewers, and documentation were written jointly

---

### Status

This repository is a **work in progress**. Parsers are being actively developed and validated. Documentation is being added progressively.

**Current focus:** Effelsberg/MBFITS parser validation with real Holmberg 1 HI data.

---

### Contact

Francesco Di Giovanni
HALO Project — Bolzano/Bozen, Italy
halo.observatory.bz [at] gmail [dot] com

If you notice any errors or have suggestions for additions, please feel free to open an issue or contact us directly.

*"Sky was just the beginning."*
