#!/usr/bin/env python3
"""
plot_feasts.py — Full visualization suite for FEASTS FAST HI data cubes.

Produces all scientifically meaningful plots from a FEASTS reduced cube:
  1. Moment-0 map (integrated HI flux)
  2. Peak flux map
  3. Moment-1 map (intensity-weighted velocity field)
  4. Moment-2 map (velocity dispersion)
  5. NHI column density map
  6. Moment-0 + contours overlaid on Peak and Moment-1
  7. Integrated spectrum over full disk (double-horn profile)
  8. Position-velocity diagram along major axis

Usage: python3 plot_feasts.py <path_to_cube.fits>
       or: python3 plot_feasts.py  (prompts for path)

Author: Francesco Di Giovanni (Claude AI assisted) - FRO/HALO Project
"""

import sys
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from astropy.io import fits
from astropy.wcs import WCS

# NHI conversion factor [cm^-2 per (Jy/beam * km/s)]
# NHI = factor * integral(S dv), assuming optically thin
# For FAST beam: factor depends on beam solid angle
# Standard single-dish approximation: NHI [cm^-2] = 1.823e18 * integral [K km/s]
# FEASTS BUNIT = Jy/beam; conversion Jy/beam -> K requires beam area
# Beam FWHM = 3.24 arcmin -> beam solid angle -> Jy/beam to K factor
# K = Jy/beam * lambda^2 / (2 * k_B * Omega_beam)
# For FAST at 1420 MHz: ~1 Jy/beam ~ 6.1 K (approximate, beam-dependent)
JY_BEAM_TO_K = 6.1   # approximate for FAST 3.24' beam at 1420 MHz
NHI_FACTOR   = 1.823e18  # cm^-2 / (K km/s)


def load_cube(cube_path):
    with fits.open(cube_path) as hdul:
        header = hdul[0].header
        data   = hdul[0].data.astype(np.float64)  # (NAXIS3, NAXIS2, NAXIS1)
    return header, data


def build_vel_axis(header):
    naxis3 = header['NAXIS3']
    crpix3 = header['CRPIX3']
    crval3 = header['CRVAL3']
    cdelt3 = header['CDELT3']
    v = crval3 + (np.arange(1, naxis3 + 1) - crpix3) * cdelt3
    return v / 1000.0  # m/s -> km/s


def compute_moments(data, vel_kms):
    dv = abs(vel_kms[1] - vel_kms[0])
    v3d = vel_kms[:, np.newaxis, np.newaxis]

    # Positive flux only for moment calculations
    data_pos = np.where(data > 0, data, np.nan)

    mom0 = np.nansum(data_pos, axis=0) * dv                          # Jy/beam * km/s
    peak = np.nanmax(data, axis=0)                                    # Jy/beam
    flux_sum = np.nansum(data_pos, axis=0)
    mom1 = np.nansum(data_pos * v3d, axis=0) / np.where(flux_sum > 0, flux_sum, np.nan)
    mom2 = np.sqrt(
        np.nansum(data_pos * (v3d - mom1[np.newaxis, :, :]) ** 2, axis=0)
        / np.where(flux_sum > 0, flux_sum, np.nan)
    )
    nhi  = mom0 * JY_BEAM_TO_K * NHI_FACTOR                          # cm^-2

    return mom0, peak, mom1, mom2, nhi


def make_wcs2d(header):
    return WCS(header).celestial


def plot_maps(mom0, peak, mom1, mom2, nhi, wcs2d, galaxy, plots_dir):
    """Five individual maps + two overlay maps."""

    # --- Contour levels for moment-0 overlay ---
    mom0_finite = mom0[np.isfinite(mom0)]
    vmax_m0 = np.nanpercentile(mom0_finite, 99)
    contour_levels = vmax_m0 * np.array([0.05, 0.15, 0.30, 0.50, 0.70, 0.90])

    configs = [
        ('mom0',  mom0,  'Moment-0 (Jy/beam km/s)', 'inferno',  None,   None),
        ('peak',  peak,  'Peak Flux (Jy/beam)',      'hot',      None,   None),
        ('mom1',  mom1,  'Moment-1 (km/s)',          'coolwarm', None,   None),
        ('mom2',  mom2,  'Moment-2 (km/s)',          'plasma',   None,   None),
        ('nhi',   nhi,   'NHI (cm⁻²)',               'viridis',  None,   None),
    ]

    for tag, mapdata, title, cmap, vmin, vmax in configs:
        fig, ax = plt.subplots(1, 1, figsize=(7, 6),
                               subplot_kw={'projection': wcs2d})
        fig.patch.set_facecolor('black')
        ax.set_facecolor('black')
        im = ax.imshow(mapdata, origin='lower', cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation='nearest')
        cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cb.ax.yaxis.set_tick_params(color='white')
        plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
        ax.set_title(f'{galaxy} — {title}', color='white', fontsize=11)
        ax.coords[0].set_axislabel('RA (J2000)', color='white')
        ax.coords[1].set_axislabel('Dec (J2000)', color='white')
        ax.coords[0].set_ticklabel(color='white')
        ax.coords[1].set_ticklabel(color='white')
        for spine in ax.spines.values():
            spine.set_edgecolor('white')
        plt.tight_layout()
        outpath = os.path.join(plots_dir, f'{galaxy}_map_{tag}.png')
        plt.savefig(outpath, dpi=150, facecolor='black')
        plt.close()
        print(f"  Saved: {outpath}")

    # --- Overlay: Peak + Moment-0 contours ---
    fig, ax = plt.subplots(1, 1, figsize=(7, 6),
                           subplot_kw={'projection': wcs2d})
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    im = ax.imshow(peak, origin='lower', cmap='hot', interpolation='nearest')
    ax.contour(mom0, levels=contour_levels, colors='cyan', linewidths=0.7, alpha=0.8)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax.set_title(f'{galaxy} — Peak Flux + Moment-0 contours', color='white', fontsize=11)
    ax.coords[0].set_axislabel('RA (J2000)', color='white')
    ax.coords[1].set_axislabel('Dec (J2000)', color='white')
    ax.coords[0].set_ticklabel(color='white')
    ax.coords[1].set_ticklabel(color='white')
    plt.tight_layout()
    outpath = os.path.join(plots_dir, f'{galaxy}_map_peak_contours.png')
    plt.savefig(outpath, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {outpath}")

    # --- Overlay: Moment-1 + Moment-0 contours ---
    fig, ax = plt.subplots(1, 1, figsize=(7, 6),
                           subplot_kw={'projection': wcs2d})
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    im = ax.imshow(mom1, origin='lower', cmap='coolwarm', interpolation='nearest')
    ax.contour(mom0, levels=contour_levels, colors='white', linewidths=0.7, alpha=0.8)
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax.set_title(f'{galaxy} — Moment-1 + Moment-0 contours', color='white', fontsize=11)
    ax.coords[0].set_axislabel('RA (J2000)', color='white')
    ax.coords[1].set_axislabel('Dec (J2000)', color='white')
    ax.coords[0].set_ticklabel(color='white')
    ax.coords[1].set_ticklabel(color='white')
    plt.tight_layout()
    outpath = os.path.join(plots_dir, f'{galaxy}_map_mom1_contours.png')
    plt.savefig(outpath, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {outpath}")


def plot_integrated_spectrum(data, vel_kms, galaxy, plots_dir):
    """Integrated spectrum over all spatial pixels — shows double-horn profile."""
    # Sum over spatial axes, NaN-safe
    integrated = np.nansum(data, axis=(1, 2))  # Jy/beam * npix

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    ax.plot(vel_kms, integrated, color='cyan', linewidth=1.0)
    ax.axhline(0, color='white', linewidth=0.5, linestyle='--')
    ax.set_xlabel('LSR velocity (km/s)', color='white')
    ax.set_ylabel('Integrated flux (Jy/beam × pixels)', color='white')
    ax.set_title(f'{galaxy} — FEASTS FAST HI integrated spectrum', color='white')
    ax.tick_params(colors='white')
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')

    # Annotate peak
    idx = np.argmax(integrated)
    ax.annotate(
        f'Peak: {integrated[idx]:.1f}\nv={vel_kms[idx]:.1f} km/s\n{1420.405752 * (1 - vel_kms[idx] / 299792.458):.4f} MHz',
        xy=(vel_kms[idx], integrated[idx]),
        xytext=(vel_kms[idx] + 60, integrated[idx] * 0.85),
        arrowprops=dict(arrowstyle='->', color='yellow'),
        color='yellow', fontsize=9
    )
    # Secondary frequency axis (GHz)
    f0_ghz = 1.420405752  # HI rest frequency GHz
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    v_ticks = ax.get_xticks()
    f_ticks = f0_ghz * (1 - v_ticks / 299792.458)
    ax2.set_xticks(v_ticks)
    ax2.set_xticklabels([f'{f:.4f}' for f in f_ticks], fontsize=7, color='white')
    ax2.set_xlabel('Frequency (GHz)', color='white')
    ax2.tick_params(colors='white')
    ax2.spines['top'].set_color('white')
    # Secondary frequency axis (GHz)
    f0_ghz = 1.420405752  # HI rest frequency GHz
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    v_ticks = ax.get_xticks()
    f_ticks = f0_ghz * (1 - v_ticks / 299792.458)
    ax2.set_xticks(v_ticks)
    ax2.set_xticklabels([f'{f:.4f}' for f in f_ticks], fontsize=7, color='white')
    ax2.set_xlabel('Frequency (GHz)', color='white')
    ax2.tick_params(colors='white')
    ax2.spines['top'].set_color('white')
    plt.tight_layout()
    outpath = os.path.join(plots_dir, f'{galaxy}_integrated_spectrum.svg')
    plt.savefig(outpath, facecolor="black")
    plt.close()
    print(f"  Saved: {outpath}")


def plot_pv_diagram(data, vel_kms, header, galaxy, plots_dir):
    """Position-velocity diagram along RA axis through Dec centre."""
    naxis2 = header['NAXIS2']
    mid_row = naxis2 // 2
    pv = data[:, mid_row, :]  # (NAXIS3, NAXIS1)

    naxis1 = header['NAXIS1']
    crpix1 = header['CRPIX1']
    cdelt1 = header['CDELT1']
    # Offset in arcmin from centre
    pix_offset = (np.arange(1, naxis1 + 1) - crpix1) * cdelt1 * 60.0  # arcmin

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.patch.set_facecolor('black')
    ax.set_facecolor('black')
    import matplotlib.colors as mcolors
    noise_region = pv[:200, :]
    sigma = np.nanstd(noise_region)
    pv_masked = np.where(pv > 3*sigma, pv, np.nan)
    pv_pos = pv_masked[np.isfinite(pv_masked)]
    vmin = np.percentile(pv_pos, 1)
    vmax = np.percentile(pv_pos, 99.5)
    im = ax.imshow(
        pv_masked, origin='lower', aspect='auto', cmap='inferno',
        norm=mcolors.PowerNorm(gamma=0.45, vmin=vmin, vmax=vmax),
        extent=[pix_offset[0], pix_offset[-1], vel_kms[0], vel_kms[-1]]
    )
    cb = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('Jy/beam', color='white')
    cb.ax.yaxis.set_tick_params(color='white')
    plt.setp(cb.ax.yaxis.get_ticklabels(), color='white')
    ax.set_xlabel('RA offset (arcmin)', color='white')
    ax.set_ylabel('LSR velocity (km/s)', color='white')
    ax.set_title(f'{galaxy} — FEASTS FAST HI Position-Velocity (along RA)', color='white')
    ax.tick_params(colors='white')
    for spine in ax.spines.values():
        spine.set_color('white')
    plt.tight_layout()
    outpath = os.path.join(plots_dir, f'{galaxy}_pv_diagram.png')
    plt.savefig(outpath, dpi=150, facecolor='black')
    plt.close()
    print(f"  Saved: {outpath}")



def plot_centre_spectrum(data, vel_kms, header, galaxy, plots_dir):
    """Spectrum extracted from the central pixel of the galaxy."""
    naxis1 = header['NAXIS1']
    naxis2 = header['NAXIS2']
    cx = naxis1 // 2
    cy = naxis2 // 2
    spec = data[:, cy, cx]

    # Literature systemic velocity (heliocentric) and measured peak
    v_sys_lit = 657.0   # km/s heliocentric (NGC 628)
    idx_peak  = np.argmax(spec)
    v_peak    = vel_kms[idx_peak]
    s_peak    = spec[idx_peak]

    fig, ax = plt.subplots(figsize=(10, 5))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.plot(vel_kms, spec, color='steelblue', linewidth=1.0)
    ax.axvline(v_sys_lit, color='gray',  linewidth=1.2, linestyle='--',
               label=f'v_sys={v_sys_lit:.0f} km/s (literature, heliocentric)')
    ax.axvline(v_peak,    color='red',   linewidth=1.2, linestyle='--',
               label=f'Peak measured={v_peak:.1f} km/s (barycentric)')
    ax.annotate(
        f'{s_peak:.4f} Jy/beam\n{1420.405752 * (1 - v_peak / 299792.458):.4f} MHz',
        xy=(v_peak, s_peak),
        xytext=(v_peak + 60, s_peak * 0.85),
        arrowprops=dict(arrowstyle='->', color='red'),
        color='red', fontsize=9
    )
    ax.set_xlabel('LSR velocity (km/s)')
    ax.set_ylabel('Jy/beam')
    ax.set_title(f'{galaxy} - FEASTS FAST HI spectrum (centre)')
    ax.legend(fontsize=8)
    # Secondary frequency axis (GHz)
    f0_ghz = 1.420405752  # HI rest frequency GHz
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    v_ticks = ax.get_xticks()
    f_ticks = f0_ghz * (1 - v_ticks / 299792.458)
    ax2.set_xticks(v_ticks)
    ax2.set_xticklabels([f'{f:.4f}' for f in f_ticks], fontsize=7)
    ax2.set_xlabel('Frequency (GHz)')
    plt.tight_layout()
    outpath = os.path.join(plots_dir, f'{galaxy}_centre_spectrum.svg')
    plt.savefig(outpath, facecolor='white')
    plt.close()
    print(f"  Saved: {outpath}")

def main():
    if len(sys.argv) > 1:
        cube_path = sys.argv[1]
    else:
        cube_path = input("Path to FEASTS cube FITS file: ").strip()
    cube_path = os.path.expanduser(cube_path)
    if not os.path.isfile(cube_path):
        print(f"ERROR: file not found: {cube_path}")
        sys.exit(1)

    print(f"Loading {cube_path} ...")
    header, data = load_cube(cube_path)

    galaxy    = header.get('OBJECT', os.path.basename(cube_path).split('-')[0])
    vel_kms   = build_vel_axis(header)
    wcs2d     = make_wcs2d(header)

    # Output directory auto-derived from script location
    surveys_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'Surveys')
    plots_dir   = os.path.join(surveys_dir, galaxy, 'Plots')
    os.makedirs(plots_dir, exist_ok=True)

    print("Computing moments ...")
    mom0, peak, mom1, mom2, nhi = compute_moments(data, vel_kms)

    print("Generating maps ...")
    plot_maps(mom0, peak, mom1, mom2, nhi, wcs2d, galaxy, plots_dir)

    print("Generating integrated spectrum ...")
    plot_integrated_spectrum(data, vel_kms, galaxy, plots_dir)
    print("Generating centre spectrum ...")
    plot_centre_spectrum(data, vel_kms, header, galaxy, plots_dir)

    print("Generating PV diagram ...")
    plot_pv_diagram(data, vel_kms, header, galaxy, plots_dir)

    print(f"\nAll plots saved in: {plots_dir}")

    # Open all plots
    for f in sorted(os.listdir(plots_dir)):
        if f.endswith('.png') or f.endswith('.svg'):
            subprocess.Popen(['xdg-open', os.path.join(plots_dir, f)])


if __name__ == '__main__':
    main()