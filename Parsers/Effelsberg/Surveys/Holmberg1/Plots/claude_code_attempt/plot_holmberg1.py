#!/usr/bin/env python3
# plot_holmberg1.py
# Plotting script for Effelsberg Holmberg I .fro files (FRO HDF5 v17)
# Mean spectra split by DATAPAR PHASE (flag_cal) and position-switched
# ON/OFF (flag_track), pol1/pol2 kept as separate lines.
# Produces: (1) full band overview (frequency, downsampled for readability)
#           (2) HI zoom (velocity, 1400-1440 MHz, full resolution)
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy
# 2026-08-10

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py
import subprocess

RESTFREQ_HZ = 1420.405750e6  # HI rest frequency
FULLBAND_DOWNSAMPLE = 64     # bin factor for the full-band overview only
VLSR_TARGET_KMS = 150.0      # Holmberg I systemic velocity (LSR)
VLSR_HALF_WINDOW_KMS = 500.0 # HI zoom half-window around VLSR_TARGET_KMS
YLIM_PERCENTILES = (1, 99)   # y-axis clip to suppress RFI spikes in the view
RFI_EXCLUDE_MHZ = (1450.0, 1470.0)  # broadband RFI, excluded from the
                                     # full-band y-limit calculation only

_base = os.path.expanduser('~/FRO/Parsers/Effelsberg/Surveys/Holmberg1')
FRO_PATH = os.path.join(_base, 'Holmberg1_FRO',
                         'EFFBG_2026-06-23_6723_02-26_P217mm-EDD_bb1.fro')
PLOTS_DIR = os.path.join(_base, 'Plots')

GROUPS = [(1, 1), (1, 0), (2, 1), (2, 0)]  # (phase, flag_track)
GROUP_LABELS = {(1, 1): 'Phase 1 — ON', (1, 0): 'Phase 1 — OFF',
                (2, 1): 'Phase 2 — ON', (2, 0): 'Phase 2 — OFF'}
GROUP_COLORS = {(1, 1): '#1f77b4', (1, 0): '#17becf',
                (2, 1): '#d62728', (2, 0): '#ff7f0e'}
POL_STYLE = {'pol1': '-', 'pol2': '--'}


def safe_path(directory, fname):
    """Return a path that does not overwrite an existing file.
    If directory/fname exists, appends _1, _2, ... before the extension."""
    base, ext = os.path.splitext(fname)
    candidate = os.path.join(directory, fname)
    counter = 1
    while os.path.exists(candidate):
        candidate = os.path.join(directory, f"{base}_{counter}{ext}")
        counter += 1
    return candidate


def vel_axis(freq_hz):
    return (RESTFREQ_HZ - freq_hz) / RESTFREQ_HZ * 2.998e5


def downsample_mean(freq_hz, spec, factor):
    n_trim = (freq_hz.size // factor) * factor
    freq_ds = freq_hz[:n_trim].reshape(-1, factor).mean(axis=1)
    spec_ds = spec[:n_trim].reshape(-1, factor).mean(axis=1)
    return freq_ds, spec_ds


def compute_group_means(fro_path, block_size=300):
    """Streams raw_spectrum in row blocks and accumulates per-(phase,
    flag_track, polarization) sums, avoiding loading the full ~8 GB
    dataset into memory at once."""
    with h5py.File(fro_path, 'r') as f:
        freq_hz = f['Spectra/frequency_axis_hz'][:]
        pol = f['Spectra/polarization'][:].astype(str)
        track = f['Quality/flag_track'][:]
        cal = np.round(f['Quality/flag_cal'][:]).astype(int)
        dset = f['Spectra/raw_spectrum']
        n, channels = dset.shape

        pols = sorted(set(pol))
        sums = {(phase, on, p): np.zeros(channels, dtype=np.float64)
                for phase, on in GROUPS for p in pols}
        counts = {key: 0 for key in sums}

        for start in range(0, n, block_size):
            end = min(start + block_size, n)
            block = dset[start:end, :]
            block_cal = cal[start:end]
            block_track = track[start:end]
            block_pol = pol[start:end]
            for phase, on in GROUPS:
                for p in pols:
                    mask = (block_cal == phase) & (block_track == on) & (block_pol == p)
                    n_hit = int(mask.sum())
                    if n_hit:
                        sums[(phase, on, p)] += block[mask].sum(axis=0, dtype=np.float64)
                        counts[(phase, on, p)] += n_hit

        attrs = {
            'facility': f['Source'].attrs.get('facility', 'Effelsberg'),
            'instrument': f['Source'].attrs.get('instrument', ''),
            'dataset': f['Source'].attrs.get('dataset', ''),
            'target': f['Observation'].attrs.get('target', ''),
            'calnote': f['Observation'].attrs.get('temperature_calibration_note', ''),
        }

    means = {key: (sums[key] / counts[key] if counts[key] > 0
                   else np.full(channels, np.nan))
              for key in sums}
    return freq_hz, means, counts, pols, attrs


def set_percentile_ylim(ax, arrays):
    """Clips the y-axis view to the 1st/99th percentile of the plotted
    data so RFI spikes don't drown out the HI signal. Data values
    themselves are untouched, only the axis range."""
    stacked = np.concatenate([a[np.isfinite(a)] for a in arrays])
    if stacked.size == 0:
        return
    lo, hi = np.nanpercentile(stacked, YLIM_PERCENTILES)
    margin = 0.1 * (hi - lo) if hi > lo else 1.0
    ax.set_ylim(lo - margin, hi + margin)


def set_fullband_ylim(ax, arrays):
    """Full-band y-limit: lower bound from the 1st percentile (RFI band
    already excluded from `arrays`), upper bound from the actual max of
    that same RFI-excluded data, so no non-RFI spike is cut off at the
    top border."""
    stacked = np.concatenate([a[np.isfinite(a)] for a in arrays])
    if stacked.size == 0:
        return
    lo = np.nanpercentile(stacked, YLIM_PERCENTILES[0])
    hi = np.nanmax(stacked)
    margin = 0.05 * (hi - lo) if hi > lo else 1.0
    ax.set_ylim(lo - margin, hi + margin)


def plot_holmberg1(fro_path, out_dir):
    print('Computing PHASE/ON-OFF mean spectra (streamed)...')
    freq_hz, means, counts, pols, attrs = compute_group_means(fro_path)
    tag = os.path.splitext(os.path.basename(fro_path))[0]
    title = f"Effelsberg 100m — {attrs['target']} — {attrs['instrument']}"
    plots = []

    # Figure 1: full band overview (downsampled), frequency axis
    fig, ax = plt.subplots(figsize=(16, 8), facecolor='black')
    ax.set_facecolor('black')
    fullband_arrays = []
    rfi_mask_ds = None
    for phase, on in GROUPS:
        for p in pols:
            key = (phase, on, p)
            freq_ds, spec_ds = downsample_mean(freq_hz, means[key], FULLBAND_DOWNSAMPLE)
            if rfi_mask_ds is None:
                freq_ds_mhz = freq_ds / 1e6
                rfi_mask_ds = (freq_ds_mhz < RFI_EXCLUDE_MHZ[0]) | (freq_ds_mhz > RFI_EXCLUDE_MHZ[1])
            ax.plot(freq_ds / 1e6, spec_ds, color=GROUP_COLORS[(phase, on)],
                    linestyle=POL_STYLE[p], linewidth=0.7,
                    label=f'{GROUP_LABELS[(phase, on)]} ({p}, n={counts[key]})')
            fullband_arrays.append(spec_ds[rfi_mask_ds])
    ax.axvline(RESTFREQ_HZ / 1e6, color='yellow', linewidth=0.8, linestyle='--',
               alpha=0.7, label=f'HI {RESTFREQ_HZ/1e6:.3f} MHz')
    set_fullband_ylim(ax, fullband_arrays)
    ax.set_xlabel('Frequency (MHz)', color='white', fontsize=11)
    ax.set_ylabel('T$_A$ (K)', color='white', fontsize=11)
    ax.set_title(f'{title} — full band, PHASE/ON-OFF', color='white', fontsize=13, pad=10)
    ax.tick_params(colors='white')
    ax.spines[:].set_edgecolor('white')
    ax.legend(fontsize=7, facecolor='#222', edgecolor='white', labelcolor='white', ncol=2)
    ax.grid(True, color='#333', linewidth=0.4)
    fig.text(0.5, -0.02, attrs['calnote'], ha='center', fontsize=7, color='gray')
    fig.tight_layout()
    out1 = safe_path(out_dir, f'{tag}_fullband_PHASE_ONOFF.svg')
    fig.savefig(out1, bbox_inches='tight', facecolor='black', edgecolor='none')
    plt.close(fig)
    print(f'  Written: {out1}')
    plots.append(out1)

    # Figure 2: HI zoom around Holmberg I systemic velocity, full resolution
    vel_full = vel_axis(freq_hz)
    vmin = VLSR_TARGET_KMS - VLSR_HALF_WINDOW_KMS
    vmax = VLSR_TARGET_KMS + VLSR_HALF_WINDOW_KMS
    hi_mask = (vel_full >= vmin) & (vel_full <= vmax)
    if np.sum(hi_mask) > 10:
        vel = vel_full[hi_mask]
        fig, ax = plt.subplots(figsize=(16, 8), facecolor='black')
        ax.set_facecolor('black')
        zoom_arrays = []
        for phase, on in GROUPS:
            for p in pols:
                key = (phase, on, p)
                spec_zoom = means[key][hi_mask]
                ax.plot(vel, spec_zoom, color=GROUP_COLORS[(phase, on)],
                        linestyle=POL_STYLE[p], linewidth=0.8,
                        label=f'{GROUP_LABELS[(phase, on)]} ({p}, n={counts[key]})')
                zoom_arrays.append(spec_zoom)
        ax.axvline(0, color='yellow', linewidth=0.8, linestyle='--',
                   alpha=0.7, label='v=0 km/s')
        ax.axvline(VLSR_TARGET_KMS, color='magenta', linewidth=0.8, linestyle=':',
                   alpha=0.7, label=f'HOI systemic v={VLSR_TARGET_KMS:.0f} km/s')
        set_percentile_ylim(ax, zoom_arrays)
        ax.set_xlabel('LSR Velocity (km/s)', color='white', fontsize=11)
        ax.set_ylabel('T$_A$ (K)', color='white', fontsize=11)
        ax.set_title(f'{title} — HI zoom (v={VLSR_TARGET_KMS:.0f}±{VLSR_HALF_WINDOW_KMS:.0f} km/s), PHASE/ON-OFF',
                     color='white', fontsize=13, pad=10)
        ax.tick_params(colors='white')
        ax.spines[:].set_edgecolor('white')
        ax.legend(fontsize=7, facecolor='#222', edgecolor='white', labelcolor='white', ncol=2)
        ax.grid(True, color='#333', linewidth=0.4)
        fig.text(0.5, -0.02, attrs['calnote'], ha='center', fontsize=7, color='gray')
        fig.tight_layout()
        out2 = safe_path(out_dir, f'{tag}_HI_zoom_PHASE_ONOFF.svg')
        fig.savefig(out2, bbox_inches='tight', facecolor='black', edgecolor='none')
        plt.close(fig)
        print(f'  Written: {out2}')
        plots.append(out2)

    return plots


def main():
    fro_path = os.path.expanduser(sys.argv[1]) if len(sys.argv) > 1 else FRO_PATH
    out_dir = os.path.expanduser(sys.argv[2]) if len(sys.argv) > 2 else PLOTS_DIR
    os.makedirs(out_dir, exist_ok=True)

    print(f'Plotting: {os.path.basename(fro_path)}')
    plots = plot_holmberg1(fro_path, out_dir)

    for p in plots:
        subprocess.Popen(['xdg-open', p],
                         stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print('\nDone.')


if __name__ == '__main__':
    main()
