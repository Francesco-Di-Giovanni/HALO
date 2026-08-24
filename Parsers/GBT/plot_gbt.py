#!/usr/bin/env python3
# plot_gbt.py
# Plotting script for GBT SDFITS .fro files (FRO HDF5 v17)
# Produces: averaged spectrum (PNG+SVG), all spectra overlay (PNG), interactive HTML
# HALO Project — Francesco Di Giovanni, Bolzano-Bozen, Italy
# 2026-08-23

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import h5py


def vel_axis(freq_hz, restfreq_hz):
    return (restfreq_hz - freq_hz) / restfreq_hz * 2.998e5


def load_fro(path):
    with h5py.File(path, 'r') as f:
        freq_hz   = f['Spectra/frequency_axis_hz'][:]
        spectra   = f['Spectra/raw_spectrum'][:]
        pols      = f['Spectra/polarization'][:].astype(str)
        obs       = f['Observation'].attrs
        target    = str(obs.get('target', ''))
        centre_hz = float(obs.get('center_frequency_hz', 0.0))
        obs_id    = f['Source'].attrs.get('observation_id', '')
        if not target:
            parts = str(obs_id).split('_')
            target = parts[3] if len(parts) > 3 else os.path.basename(path)
    return freq_hz, spectra, pols, target, centre_hz


def band_label(centre_hz):
    if centre_hz < 5e9:
        return 'HI 21cm'
    else:
        return f'H2O 22 GHz'


def plot_fro(fro_path, out_dir, plotly_js):
    freq_hz, spectra, pols, target, centre_hz = load_fro(fro_path)
    restfreq_hz = centre_hz if centre_hz > 0 else 1420.405750e6
    vel = vel_axis(freq_hz, restfreq_hz)
    tag = os.path.splitext(os.path.basename(fro_path))[0]
    unique_pols = list(dict.fromkeys(pols))
    colors = {'XX': '#1f77b4', 'YY': '#d62728',
              'RR': '#1f77b4', 'LL': '#d62728', 'I': '#2ca02c'}
    title = f'GBT — {target} — {band_label(centre_hz)}'

    # ── Figure 1: averaged spectrum (PNG + SVG) ───────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 6), facecolor='white')
    ax.set_facecolor('white')
    avg_by_pol = {}
    for pol in unique_pols:
        mask = pols == pol
        avg = np.nanmean(spectra[mask], axis=0)
        avg_by_pol[pol] = avg
        ax.plot(vel, avg, color=colors.get(pol, '#7f7f7f'),
                linewidth=0.8, label=pol)
    ax.axvline(0, color='gray', linewidth=0.7, linestyle='--',
               alpha=0.6, label='v=0 km/s')
    ax.set_xlabel('LSR Velocity (km/s)', color='black', fontsize=9)
    ax.set_ylabel('T$_A$ (K)', color='black', fontsize=9)
    ax.set_title(title, color='black', fontsize=11, pad=8)
    ax.tick_params(colors='black')
    ax.spines[:].set_edgecolor('black')
    ax.legend(fontsize=8, facecolor='white', edgecolor='black', labelcolor='black')
    ax.grid(True, color='#cccccc', linewidth=0.4)
    # Second x-axis: frequency in MHz
    def vel_to_freq(v): return restfreq_hz/1e6 * (1 - v/2.998e5)
    def freq_to_vel(f): return (1 - f/(restfreq_hz/1e6)) * 2.998e5
    ax2 = ax.secondary_xaxis(-0.22, functions=(vel_to_freq, freq_to_vel))
    ax2.set_xlabel('Frequency (MHz)', color='black', fontsize=9)
    ax2.tick_params(colors='black')
    fig.subplots_adjust(bottom=0.22)
    fig.tight_layout(rect=[0, 0.18, 1, 1])
    out_png = os.path.join(out_dir, f'{tag}_avg_spectrum.png')
    out_svg = os.path.join(out_dir, f'{tag}_avg_spectrum.svg')
    fig.savefig(out_png, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    fig.savefig(out_svg, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f'  Written: {out_png}')
    print(f'  Written: {out_svg}')

    # ── Figure 2: all spectra overlay (PNG) ───────────────────────────────────
    fig, ax = plt.subplots(figsize=(16, 6), facecolor='white')
    ax.set_facecolor('white')
    for spec, pol in zip(spectra, pols):
        ax.plot(vel, spec, color=colors.get(pol, '#7f7f7f'),
                linewidth=0.4, alpha=0.6)
    ax.axvline(0, color='gray', linewidth=0.7, linestyle='--', alpha=0.6)
    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=colors.get(p, '#7f7f7f'),
                      linewidth=1.5, label=p) for p in unique_pols]
    handles.append(Line2D([0], [0], color='gray', linewidth=1,
                          linestyle='--', label='v=0 km/s'))
    ax.legend(handles=handles, fontsize=9, facecolor='white',
              edgecolor='black', labelcolor='black')
    ax.set_xlabel('LSR Velocity (km/s)', color='black', fontsize=9)
    ax.set_ylabel('T$_A$ (K)', color='black', fontsize=9)
    ax.set_title(f'{title} — all spectra', color='black', fontsize=11, pad=8)
    ax.tick_params(colors='black')
    ax.spines[:].set_edgecolor('black')
    ax.grid(True, color='#cccccc', linewidth=0.4)
    ax2 = ax.secondary_xaxis(-0.22, functions=(vel_to_freq, freq_to_vel))
    ax2.set_xlabel('Frequency (MHz)', color='black', fontsize=9)
    ax2.tick_params(colors='black')
    fig.subplots_adjust(bottom=0.22)
    fig.tight_layout(rect=[0, 0.18, 1, 1])
    out_all = os.path.join(out_dir, f'{tag}_all_spectra.png')
    fig.savefig(out_all, dpi=150, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close(fig)
    print(f'  Written: {out_all}')

    # ── Figure 3: interactive HTML ────────────────────────────────────────────
    out_html = None
    try:
        import plotly.graph_objects as go
        fig_html = go.Figure()
        color_map = {'XX': '#1f77b4', 'YY': '#d62728',
                     'RR': '#1f77b4', 'LL': '#d62728'}
        first_pol = unique_pols[0]
        xs_js = list(vel.tolist())
        ys_js = list(avg_by_pol[first_pol].tolist())
        ys2_js = list(avg_by_pol[unique_pols[1]].tolist()) if len(unique_pols) > 1 else ys_js
        pol1_name = unique_pols[0]
        pol2_name = unique_pols[1] if len(unique_pols) > 1 else unique_pols[0]
        for pol in unique_pols:
            fig_html.add_trace(go.Scatter(
                x=vel, y=avg_by_pol[pol],
                mode='lines', name=pol,
                line=dict(color=color_map.get(pol, '#7f7f7f'), width=0.8),
                hoverinfo='none',
            ))
        fig_html.add_vline(x=0, line=dict(color='gray', width=1,
                           dash='dash'), opacity=0.6)
        fig_html.update_layout(
            title=dict(text=title, font=dict(color='black', size=14)),
            xaxis=dict(
                title='LSR Velocity (km/s)',
                color='black', gridcolor='rgba(180,180,180,0.4)',
                minor=dict(ticks='inside', showgrid=False),
            ),
            yaxis=dict(
                title='T<sub>A</sub> (K)',
                color='black', gridcolor='rgba(180,180,180,0.4)',
                minor=dict(ticks='inside', showgrid=False),
            ),
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(color='black'),
            hovermode='closest',
            dragmode='zoom',
            autosize=True,
            hoverlabel=dict(namelength=0),
            legend=dict(bgcolor='rgba(240,240,240,0.8)', bordercolor='black',
                        borderwidth=1, font=dict(color='black')),
            annotations=[dict(
                x=0.01, y=0.99,
                xref='paper', yref='paper',
                xanchor='left', yanchor='top',
                text='v: — km/s<br>freq: — MHz<br>pol1: —<br>pol2: —',
                showarrow=False,
                bgcolor='rgba(240,240,240,0.85)',
                bordercolor='gray',
                borderwidth=1,
                font=dict(size=12, family='monospace', color='black'),
                name='info_box',
            )],
        )
        fs_js = list((restfreq_hz/1e6 * (1 - vel/2.998e5)).tolist())
        info_box_js = (
            "<script>"
            "var xs=" + str(xs_js) + ";"
            "var ys=" + str(ys_js) + ";"
            "var ys2=" + str(ys2_js) + ";"
            "var fs=" + str(fs_js) + ";"
            "var pol1='" + pol1_name + "';"
            "var pol2='" + pol2_name + "';"
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
            "      'annotations[0].text':'<b>v</b>: '+xi.toFixed(1)+' km/s<br><b>freq</b>: '+fs[best].toFixed(3)+' MHz<br><b>'+pol1+'</b>: '+yi.toFixed(4)+' K<br><b>'+pol2+'</b>: '+ys2[best].toFixed(4)+' K',"
            "      shapes:["
            "        {type:'line',x0:xi,x1:xi,y0:0,y1:1,yref:'paper',line:{color:'gray',width:1,dash:'dot'}},"
            "        {type:'line',x0:0,x1:1,xref:'paper',y0:yi,y1:yi,line:{color:'gray',width:1,dash:'dot'}}"
            "      ]"
            "    });"
            "  });"
            "  myPlot.addEventListener('mouseleave',function(){"
            "    Plotly.relayout(myPlot,{'annotations[0].text':'v: — km/s<br>freq: — MHz<br>pol1: —<br>pol2: —',shapes:[]});"
            "  });"
            "},500);"
            "</script>"
        )
        out_html = os.path.join(out_dir, f'{tag}_avg_spectrum.html')
        if plotly_js and os.path.isfile(plotly_js):
            with open(plotly_js, 'r') as pf:
                js_src = pf.read()
            html_content = fig_html.to_html(full_html=True, include_plotlyjs=False)
            html_content = html_content.replace(
                '<head>', f'<head><script>{js_src}</script>', 1)
        else:
            html_content = fig_html.to_html(full_html=True, include_plotlyjs=True)
        html_content = html_content.replace('</body>', info_box_js + '</body>')
        with open(out_html, 'w') as f_html:
            f_html.write(html_content)
        print(f'  Written: {out_html}')
    except Exception as e:
        print(f'  HTML plot skipped: {e}')

    return [p for p in [out_png, out_svg, out_all, out_html] if p]


def main():
    plotly_js = os.path.expanduser(
        '~/.local/lib/python3.14/site-packages/plotly/package_data/plotly.min.js')
    if not os.path.isfile(plotly_js):
        plotly_js = None

    if len(sys.argv) < 2:
        base = os.path.expanduser('~/FRO/Parsers/GBT/Surveys/')
        fro_files = sorted(glob.glob(
            os.path.join(base, '*', 'GBT_FRO_*', '*.h5')))
        if not fro_files:
            print('No .fro files found.')
            sys.exit(1)
        print(f'Found {len(fro_files)} .fro file(s):')
        for i, f in enumerate(fro_files):
            print(f'  [{i}] {os.path.basename(f)}')
        choice = input('Select file number (or Enter for all): ').strip()
        selected = fro_files if choice == '' else [fro_files[int(choice)]]
    else:
        selected = [os.path.expanduser(sys.argv[1])]

    import subprocess
    all_plots = []
    for fro_path in selected:
        fro_dir = os.path.dirname(fro_path)
        obs_dir = os.path.dirname(fro_dir)
        plots_dirname = os.path.basename(fro_dir).replace('GBT_FRO_', 'GBT_Plots_')
        out_dir = os.path.expanduser(sys.argv[2]) if len(sys.argv) >= 3 else \
                  os.path.join(obs_dir, plots_dirname)
        os.makedirs(out_dir, exist_ok=True)
        print(f'\nPlotting: {os.path.basename(fro_path)}')
        all_plots.extend(plot_fro(fro_path, out_dir, plotly_js))
    for p in all_plots:
        if p and os.path.isfile(p):
            subprocess.Popen(['xdg-open', p],
                             stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    print('\nDone.')


if __name__ == '__main__':
    main()
