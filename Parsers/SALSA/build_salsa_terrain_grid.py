# build_salsa_terrain_grid.py
# Build a 2D GLon x GLat peak-temperature grid from the SALSA galactic survey
# (SALSA_Survey_Galactic_5deg_120s), for 3D terrain visualization.
#
# For each FITS pointing, computes two peak values:
#   A) masked only (DC spike + known RFI window excluded)
#   B) masked + per-spectrum baseline subtracted (removes passband slope/ripple)
#
# Output: NPZ file with glon_values, glat_values (1D axes) and peak_grid_a,
# peak_grid_b (2D arrays, NaN where no observation exists).
#
# FRO/HALO Project - Francesco Di Giovanni, Bolzano

import os
import glob
import numpy as np
from astropy.io import fits

SURVEY_DIR = os.path.expanduser(
    "~/FRO/Parsers/SALSA/Surveys/SALSA_Survey_Galactic_5deg_120s"
)
OUTPUT_NPZ = os.path.join(SURVEY_DIR, "salsa_terrain_grid.npz")

REST_FREQ_HZ = 1420.405e6
C_KMS = 299792.458

# RFI window to exclude, in km/s (LSR velocity)
RFI_VMIN_KMS = 80.0
RFI_VMAX_KMS = 95.0

# DC spike: exclude central channel +/- margin
DC_MARGIN_CHANNELS = 2

# Baseline: fraction of channels at each edge of the band used for
# baseline subtraction (avoids the HI line region in the center)
BASELINE_EDGE_FRACTION = 0.07

# Grid step, must match the survey's actual step (deg)
GRID_STEP_DEG = 5.0


def build_frequency_axis(header):
    n = int(header["NAXIS1"])
    crpix1 = header["CRPIX1"]
    crval1 = header["CRVAL1"]
    cdelt1 = header["CDELT1"]
    pixels = np.arange(1, n + 1)
    return crval1 + (pixels - crpix1) * cdelt1


def compute_peaks(spectrum, velocity):
    """Return (peak_a, peak_b) for one spectrum.
    peak_a: DC + RFI masked only.
    peak_b: DC + RFI masked, plus baseline-subtracted.
    """
    n = len(spectrum)
    spec = spectrum.astype(float).copy()

    mask = np.ones(n, dtype=bool)

    center = n // 2
    lo = max(0, center - DC_MARGIN_CHANNELS)
    hi = min(n, center + DC_MARGIN_CHANNELS + 1)
    mask[lo:hi] = False

    rfi_mask = (velocity >= RFI_VMIN_KMS) & (velocity <= RFI_VMAX_KMS)
    mask[rfi_mask] = False

    spec_a = np.where(mask, spec, np.nan)
    peak_a = float(np.nanmax(spec_a))

    edge_n = int(n * BASELINE_EDGE_FRACTION)
    baseline_idx = np.zeros(n, dtype=bool)
    baseline_idx[:edge_n] = True
    baseline_idx[-edge_n:] = True
    baseline_idx &= mask

    baseline_level = np.mean(spec[baseline_idx])
    spec_b = spec - baseline_level
    spec_b = np.where(mask, spec_b, np.nan)
    spec_b = np.clip(spec_b, 0, None)
    peak_b = float(np.nanmax(spec_b))

    return peak_a, peak_b


def main():
    fits_files = sorted(glob.glob(os.path.join(SURVEY_DIR, "*.fits")))
    print(f"Found {len(fits_files)} FITS files in survey directory.")

    records = []
    skipped = 0

    for fp in fits_files:
        try:
            with fits.open(fp) as hdul:
                header = hdul[0].header
                data = hdul[0].data
                glon = header.get("CRVAL2", np.nan)
                glat = header.get("CRVAL3", np.nan)

                if np.isnan(glon) or np.isnan(glat):
                    print(f"  SKIP (no GLon/GLat): {os.path.basename(fp)}")
                    skipped += 1
                    continue

                freq_hz = build_frequency_axis(header)
                velocity = C_KMS * (REST_FREQ_HZ - freq_hz) / REST_FREQ_HZ

                spectrum = np.squeeze(data)
                if spectrum.ndim != 1:
                    print(f"  SKIP (unexpected shape {data.shape}): {os.path.basename(fp)}")
                    skipped += 1
                    continue

                peak_a, peak_b = compute_peaks(spectrum, velocity)
                records.append((glon, glat, peak_a, peak_b))

        except Exception as e:
            print(f"  ERROR reading {os.path.basename(fp)}: {e}")
            skipped += 1

    print(f"Processed {len(records)} pointings, skipped {skipped}.")

    if not records:
        print("No valid records, aborting.")
        return

    records = np.array(records)
    glon_all = records[:, 0]
    glat_all = records[:, 1]
    peak_a_all = records[:, 2]
    peak_b_all = records[:, 3]

    glon_rounded = np.round(glon_all / GRID_STEP_DEG) * GRID_STEP_DEG
    glat_rounded = np.round(glat_all / GRID_STEP_DEG) * GRID_STEP_DEG

    glon_values = np.unique(glon_rounded)
    glat_values = np.unique(glat_rounded)

    peak_grid_a = np.full((len(glat_values), len(glon_values)), np.nan)
    peak_grid_b = np.full((len(glat_values), len(glon_values)), np.nan)

    for gl, gb, pa, pb in zip(glon_rounded, glat_rounded, peak_a_all, peak_b_all):
        gi = np.where(glon_values == gl)[0][0]
        bi = np.where(glat_values == gb)[0][0]
        if np.isnan(peak_grid_a[bi, gi]) or pa > peak_grid_a[bi, gi]:
            peak_grid_a[bi, gi] = pa
        if np.isnan(peak_grid_b[bi, gi]) or pb > peak_grid_b[bi, gi]:
            peak_grid_b[bi, gi] = pb

    n_cells = len(glon_values) * len(glat_values)
    n_filled = int(np.sum(~np.isnan(peak_grid_a)))
    print(f"Grid: {len(glon_values)} GLon x {len(glat_values)} GLat = {n_cells} cells, "
          f"{n_filled} filled ({n_cells - n_filled} gaps).")

    np.savez(
        OUTPUT_NPZ,
        glon_values=glon_values,
        glat_values=glat_values,
        peak_grid_a=peak_grid_a,
        peak_grid_b=peak_grid_b,
    )
    print(f"Saved: {OUTPUT_NPZ}")


if __name__ == "__main__":
    main()
