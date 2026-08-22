#!/usr/bin/env python3
"""
plot_things.py — THINGS survey visualizer for HALO
Produces 5 plots from THINGS data products:
  1. Integrated HI spectrum (from data cube)
  2. Moment 0 map (integrated HI flux)
  3. Moment 1 map (velocity field)
  4. Moment 2 map (velocity dispersion)
  5. Channel maps (grid of velocity slices)

Usage:
    python3 plot_things.py <galaxy_directory>
    e.g.: python3 plot_things.py ~/FRO/Parsers/THINGS/Surveys/NGC_5055/

Project: HALO (Hydrogen Atomic Line Observatory)
Author: Francesco Di Giovanni (Bolzano, Italy)
Created: 2026-08-21 (Parte 67)
"""

import sys
import os
import glob
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import PowerNorm
from pathlib import Path

import warnings
from astropy.io import fits
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
warnings.filterwarnings('ignore', category=AstropyWarning)


def center_window(fig):
    """Center a matplotlib TkAgg window on screen."""
    manager = fig.canvas.manager
    manager.canvas.draw()
    fig.canvas.flush_events()
    win = manager.window
    win.update_idletasks()
    sw = win.winfo_screenwidth()
    sh = win.winfo_screenheight()
    fw = win.winfo_width()
    fh = win.winfo_height()
    x = (sw - fw) // 2
    y = (sh - fh) // 2
    win.geometry(f'+{x}+{y}')

# ── HI rest frequency ────────────────────────────────────────────────────────
HI_REST_HZ = 1420405751.768
C_KMS = 299792.458  # km/s


def find_fits(directory, pattern):
    """Find a FITS file matching pattern in directory."""
    matches = glob.glob(str(Path(directory) / pattern))
    if not matches:
        return None
    return matches[0]


def load_cube(path):
    """Load a FITS cube, squeeze degenerate axes."""
    with fits.open(path) as hdul:
        hdr  = hdul[0].header.copy()
        data = hdul[0].data.copy()
    if data.ndim == 4:
        data = data.squeeze(axis=0)
    return hdr, data


def load_map(path):
    """Load a 2D FITS map."""
    with fits.open(path) as hdul:
        hdr  = hdul[0].header.copy()
        data = hdul[0].data.copy()
    while data.ndim > 2:
        data = data.squeeze(axis=0)
    return hdr, data


def get_velocity_axis(hdr, n_chan):
    """Return velocity axis in km/s from FITS header (FELO-HEL or VRAD)."""
    ctype = hdr.get('CTYPE3', '').upper().strip()
    crpix = float(hdr.get('CRPIX3', 1.0))
    crval = float(hdr.get('CRVAL3', 0.0))
    cdelt = float(hdr.get('CDELT3', 1.0))
    vel_ms = crval + (np.arange(n_chan) - (crpix - 1)) * cdelt
    return vel_ms / 1000.0  # m/s → km/s


def plot_spectrum(ax, vel_kms, spectrum, object_name, bunit):
    """Plot 1: integrated HI spectrum."""
    from matplotlib.ticker import AutoMinorLocator
    ax.plot(vel_kms, spectrum, color='steelblue', linewidth=1.0)
    ax.axhline(0, color='gray', linewidth=0.5, linestyle='--')
    ax.set_xlabel('LSR velocity (km/s)')
    ax.set_ylabel(f'Mean flux ({bunit})')
    ax.set_title(f'{object_name} — Integrated HI spectrum (THINGS)')
    ax.grid(True, alpha=0.3)
    ax.xaxis.set_minor_locator(AutoMinorLocator(5))
    ax.yaxis.set_minor_locator(AutoMinorLocator(5))
    ax.tick_params(which='minor', length=3)


def plot_moment(ax, hdr, data, title, cmap, label, vmin=None, vmax=None):
    """Plot a 2D moment map with WCS axes."""
    wcs = WCS(hdr, naxis=2)
    im = ax.imshow(data, origin='lower', cmap=cmap,
                   vmin=vmin, vmax=vmax,
                   interpolation='nearest')
    plt.colorbar(im, ax=ax, label=label, fraction=0.046, pad=0.04)
    ax.set_title(title)
    ax.set_xlabel('RA (J2000)')
    ax.set_ylabel('Dec (J2000)')
    # RA/Dec tick labels from WCS
    ra  = ax.coords[0]
    dec = ax.coords[1]
    ra.set_major_formatter('hh:mm:ss')
    dec.set_major_formatter('dd:mm:ss')
    ax.coords[0].display_minor_ticks(True)
    ax.coords[1].display_minor_ticks(True)


def plot_channel_maps(fig, cube_hdr, cube_data, vel_kms, object_name, n_cols=6):
    """Plot 5: channel maps grid."""
    n_chan = cube_data.shape[0]
    # Select evenly spaced channels
    n_show = min(24, n_chan)
    idx    = np.linspace(0, n_chan - 1, n_show, dtype=int)
    n_rows = int(np.ceil(n_show / n_cols))

    fig.suptitle(f'{object_name} — HI channel maps (THINGS)', fontsize=13)
    gs = gridspec.GridSpec(n_rows, n_cols, figure=fig,
                           hspace=0.05, wspace=0.05)

    vmax = np.nanpercentile(cube_data, 99.5)
    vmin = -0.1 * vmax

    for plot_i, chan_i in enumerate(idx):
        row = plot_i // n_cols
        col = plot_i % n_cols
        ax  = fig.add_subplot(gs[row, col],
)
        ax.imshow(cube_data[chan_i], origin='lower', cmap='inferno',
                  vmin=vmin, vmax=vmax, interpolation='nearest')
        ax.text(0.05, 0.92, f'{vel_kms[chan_i]:.0f} km/s',
                transform=ax.transAxes, color='white', fontsize=7,
                verticalalignment='top')
        ax.set_xticks([])
        ax.set_yticks([])


def main():
    if len(sys.argv) > 1:
        galaxy_dir = Path(sys.argv[1])
    else:
        import subprocess
        result = subprocess.run(
            ['zenity', '--file-selection', '--directory',
             '--title=Seleziona directory galassia THINGS',
             '--filename=' + str(Path.home() / 'FRO/Parsers/THINGS/Surveys/') + '/'],
            capture_output=True, text=True
        )
        path = result.stdout.strip()
        if not path:
            print("Nessuna directory selezionata.")
            sys.exit(0)
        galaxy_dir = Path(path)

    if not galaxy_dir.exists():
        print(f"Directory non trovata: {galaxy_dir}")
        sys.exit(1)

    # Find data files
    cube_path = find_fits(galaxy_dir, '*NA_CUBE*.FITS') or find_fits(galaxy_dir, '*NA_CUBE*.fits')
    mom0_path = find_fits(galaxy_dir, '*MOM0*.FITS') or find_fits(galaxy_dir, '*MOM0*.fits')
    mom1_path = find_fits(galaxy_dir, '*MOM1*.FITS') or find_fits(galaxy_dir, '*MOM1*.fits')
    mom2_path = find_fits(galaxy_dir, '*MOM2*.FITS') or find_fits(galaxy_dir, '*MOM2*.fits')

    if not cube_path:
        print("Cubo FITS non trovato nella directory.")
        sys.exit(1)

    print(f"Caricamento cubo: {cube_path}")
    cube_hdr, cube_data = load_cube(cube_path)
    n_chan = cube_data.shape[0]
    vel_kms = get_velocity_axis(cube_hdr, n_chan)
    object_name = cube_hdr.get('OBJECT', galaxy_dir.name).strip()
    bunit = cube_hdr.get('BUNIT', 'Jy/beam').strip()

    # Integrated spectrum
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        spectrum = np.nanmean(cube_data, axis=(1, 2))

    # Output directory
    plots_dir = galaxy_dir / 'Plots'
    plots_dir.mkdir(exist_ok=True)

    # ── Plot 1: Integrated spectrum ──────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    plot_spectrum(ax1, vel_kms, spectrum, object_name, bunit)
    fig1.tight_layout()
    center_window(fig1)
    out1 = plots_dir / f'{object_name}_HI_spectrum.png'
    fig1.savefig(out1, dpi=300, bbox_inches='tight')
    print(f"Salvato: {out1}")

    # ── Plot 1b: Interactive HTML spectrum (Plotly) ──────────────────────────
    try:
        import plotly.graph_objects as go
        fig_html = go.Figure()
        fig_html.add_trace(go.Scatter(
            x=list(vel_kms),
            y=list(spectrum),
            mode='lines',
            line=dict(color='steelblue', width=1.5),
            name='HI flux',
        ))
        fig_html.add_hline(y=0, line=dict(color='gray', width=0.8, dash='dash'))
        fig_html.update_layout(
            title=dict(text=f'{object_name} — Integrated HI spectrum (THINGS)', x=0.5),
            xaxis=dict(
                title='LSR velocity (km/s)',
                showgrid=True, gridcolor='rgba(128,128,128,0.2)',
                minor=dict(ticks='inside', showgrid=False),
            ),
            yaxis=dict(
                title=f'Mean flux ({bunit})',
                showgrid=True, gridcolor='rgba(128,128,128,0.2)',
                minor=dict(ticks='inside', showgrid=False),
            ),
            plot_bgcolor='white',
            paper_bgcolor='white',
            hovermode='closest',
            hoverlabel=dict(namelength=0),
            annotations=[dict(
                x=0.01, y=0.99,
                xref='paper', yref='paper',
                xanchor='left', yanchor='top',
                text='v: — km/s<br>flux: —',
                showarrow=False,
                bgcolor='rgba(255,255,255,0.85)',
                bordercolor='gray',
                borderwidth=1,
                font=dict(size=12, family='monospace'),
                name='info_box',
            )],
        )
        # Scale to µJy/beam for readable axis
        spectrum_ujy = spectrum * 1e6
        fig_html.data[0].y = list(spectrum_ujy)
        fig_html.update_traces(hoverinfo='none')
        fig_html.update_layout(
            yaxis_title='Mean flux (µJy/beam)',
            yaxis_tickformat='.1f',
        )
        _bunit_js = 'µJy/beam'
        info_box_js = (
            "<script>"
            "var xs=" + str(list(vel_kms.tolist())) + ";"
            "var ys=" + str(list(spectrum_ujy.tolist())) + ";"
            "setTimeout(function(){"
            "  var myPlot=document.getElementsByClassName('plotly-graph-div')[0];"
            "  myPlot.addEventListener('mousemove',function(evt){"
            "    var bb=myPlot.getBoundingClientRect();"
            "    var xl=myPlot._fullLayout.xaxis;"
            "    var mx=evt.clientX-bb.left-myPlot._fullLayout.margin.l;"
            "    var best=0,bd=Infinity;"
            "    for(var i=0;i<xs.length;i++){var d=Math.abs(xs[i]-xl.p2l(mx));if(d<bd){bd=d;best=i;}}"
            "    var xi=xs[best],yi=ys[best];"
            "    Plotly.relayout(myPlot,{"
            "      'annotations[0].text':'<b>v</b>: '+xi.toFixed(1)+' km/s<br><b>flux</b>: '+yi.toFixed(1)+' µJy/beam',"
            "      shapes:["
            "        {type:'line',x0:xi,x1:xi,y0:0,y1:1,yref:'paper',line:{color:'gray',width:1,dash:'dot'}},"
            "        {type:'line',x0:0,x1:1,xref:'paper',y0:yi,y1:yi,line:{color:'gray',width:1,dash:'dot'}}"
            "      ]"
            "    });"
            "  });"
            "  myPlot.addEventListener('mouseleave',function(){"
            "    Plotly.relayout(myPlot,{'annotations[0].text':'v: — km/s<br>flux: —',shapes:[]});"
            "  });"
            "},500);"
            "</script>"
        )
        out1_html = plots_dir / f'{object_name}_HI_spectrum.html'
        html_content = fig_html.to_html(full_html=True, include_plotlyjs=True)
        html_content = html_content.replace('</body>', info_box_js + '</body>')
        with open(out1_html, 'w') as f_html:
            f_html.write(html_content)
        print(f"Salvato: {out1_html}")
    except Exception as e:
        print(f"Errore plotly: {e}")

    # ── Plot 2: Moment 0 ─────────────────────────────────────────────────────
    if mom0_path:
        print(f"Caricamento moment 0: {mom0_path}")
        mom0_hdr, mom0_data = load_map(mom0_path)
        fig2, ax2 = plt.subplots(figsize=(8, 8), subplot_kw={'projection': WCS(mom0_hdr, naxis=2)}, constrained_layout=True)
        vmax = np.nanpercentile(mom0_data, 99.5)
        plot_moment(ax2, mom0_hdr, mom0_data,
                    f'{object_name} — HI Moment 0 (integrated flux)',
                    'inferno', mom0_hdr.get('BUNIT', 'Jy/beam km/s'),
                    vmin=0, vmax=vmax)
        out2 = plots_dir / f'{object_name}_moment0.png'
        fig2.savefig(out2, dpi=300, bbox_inches='tight')
        print(f"Salvato: {out2}")

    # ── Plot 3: Moment 1 ─────────────────────────────────────────────────────
    if mom1_path:
        print(f"Caricamento moment 1: {mom1_path}")
        mom1_hdr, mom1_data = load_map(mom1_path)
        # Convert m/s to km/s if needed
        if np.nanmax(np.abs(mom1_data)) > 1e4:
            mom1_data = mom1_data / 1000.0
        fig3, ax3 = plt.subplots(figsize=(8, 8), subplot_kw={'projection': WCS(mom1_hdr, naxis=2)}, constrained_layout=True)
        vmed = np.nanmedian(mom1_data)
        vstd = np.nanstd(mom1_data)
        plot_moment(ax3, mom1_hdr, mom1_data,
                    f'{object_name} — HI Moment 1 (velocity field)',
                    'RdBu_r', 'LSR velocity (km/s)',
                    vmin=vmed - 3.5*vstd, vmax=vmed + 3.5*vstd)
        out3 = plots_dir / f'{object_name}_moment1.png'
        fig3.savefig(out3, dpi=300, bbox_inches='tight')
        print(f"Salvato: {out3}")

    # ── Plot 4: Moment 2 ─────────────────────────────────────────────────────
    if mom2_path:
        print(f"Caricamento moment 2: {mom2_path}")
        mom2_hdr, mom2_data = load_map(mom2_path)
        if np.nanmax(np.abs(mom2_data)) > 1e4:
            mom2_data = mom2_data / 1000.0
        fig4, ax4 = plt.subplots(figsize=(8, 8), subplot_kw={'projection': WCS(mom2_hdr, naxis=2)}, constrained_layout=True)
        plot_moment(ax4, mom2_hdr, mom2_data,
                    f'{object_name} — HI Moment 2 (velocity dispersion)',
                    'plasma', 'velocity dispersion (km/s)',
                    vmin=-0.05*np.nanpercentile(mom2_data, 99), vmax=np.nanpercentile(mom2_data, 99))
        out4 = plots_dir / f'{object_name}_moment2.png'
        fig4.savefig(out4, dpi=300, bbox_inches='tight')
        print(f"Salvato: {out4}")

    # ── Plot 5a: Channel maps mosaico ───────────────────────────────────────
    fig5 = plt.figure(figsize=(14, 10))
    center_window(fig5)
    plot_channel_maps(fig5, cube_hdr, cube_data, vel_kms, object_name)
    fig5.subplots_adjust(top=0.96)
    out5 = plots_dir / f'{object_name}_channel_maps.png'
    fig5.savefig(out5, dpi=150, bbox_inches='tight')
    print(f"Salvato: {out5}")

    # ── Plot 5b: Channel maps navigabili (una finestra con frecce) ──────────
    chan_dir = plots_dir / 'Channels'
    chan_dir.mkdir(exist_ok=True)
    n_chan = cube_data.shape[0]
    n_show = min(24, n_chan)
    idx = np.linspace(0, n_chan - 1, n_show, dtype=int)
    vmax_chan = np.nanpercentile(cube_data, 99.5)
    vmin_chan = -0.1 * vmax_chan
    from matplotlib.ticker import AutoMinorLocator

    # Save individual PNGs
    for plot_i, chan_i in enumerate(idx):
        fig_c, ax_c = plt.subplots(figsize=(6, 6))
        im_c = ax_c.imshow(cube_data[chan_i], origin='lower', cmap='inferno',
                           vmin=vmin_chan, vmax=vmax_chan, interpolation='nearest')
        plt.colorbar(im_c, ax=ax_c, label=bunit, fraction=0.046, pad=0.04)
        ax_c.set_title(f'{object_name} — HI channel {vel_kms[chan_i]:.1f} km/s (THINGS)')
        ax_c.set_xlabel('pixels (RA direction)')
        ax_c.set_ylabel('pixels (Dec direction)')
        ax_c.xaxis.set_minor_locator(AutoMinorLocator())
        ax_c.yaxis.set_minor_locator(AutoMinorLocator())
        ax_c.tick_params(which='minor', length=3)
        fig_c.tight_layout()
        out_c = chan_dir / f'{object_name}_channel_{plot_i+1:02d}_{vel_kms[chan_i]:.0f}kms.png'
        fig_c.savefig(out_c, dpi=150, bbox_inches='tight')
        plt.close(fig_c)
        print(f"  Salvato: {out_c.name}")
    print(f"Channel maps singole salvate in: {chan_dir}")

    # Interactive navigator: one window, prev/next buttons
    nav_state = [0]
    fig_nav, ax_nav = plt.subplots(figsize=(7, 7))
    center_window(fig_nav)
    def draw_channel(i):
        ax_nav.clear()
        chan_i = idx[i]
        ax_nav.imshow(cube_data[chan_i], origin='lower', cmap='inferno',
                      vmin=vmin_chan, vmax=vmax_chan, interpolation='nearest')
        ax_nav.set_title(f'{object_name} — HI {vel_kms[chan_i]:.1f} km/s  [{i+1}/{n_show}]')
        ax_nav.set_xlabel('pixels (RA direction)')
        ax_nav.set_ylabel('pixels (Dec direction)')
        ax_nav.xaxis.set_minor_locator(AutoMinorLocator())
        ax_nav.yaxis.set_minor_locator(AutoMinorLocator())
        ax_nav.tick_params(which='minor', length=3)
        fig_nav.canvas.draw()

    draw_channel(0)

    print(f"\nTutti i plot salvati in: {plots_dir}")
    plt.show()


if __name__ == '__main__':
    main()
