#!/usr/bin/env python3
"""
plot_holmberg1.py - Effelsberg/Holmberg I HI interactive HTML plot script (FRO HDF5 v17)
Author: Francesco Di Giovanni (Claude AI assisted) - FRO/HALO Project
Generates 2 standalone HTML files: white background, rubber-band zoom, fixed readout.
"""

import os, sys, json, subprocess
import numpy as np
import h5py

# ── Paths ──────────────────────────────────────────────────────────────────────
FRO_PATH  = os.path.expanduser(
    '~/FRO/Parsers/Effelsberg/Surveys/Holmberg1/Holmberg1_FRO/'
    'EFFBG_2026-06-23_6723_02-26_P217mm-EDD_bb1.fro')
PLOTS_DIR = os.path.expanduser(
    '~/FRO/Parsers/Effelsberg/Surveys/Holmberg1/Plots/')

# ── Constants ──────────────────────────────────────────────────────────────────
C_KMS      = 299792.458
HI_REST_HZ = 1420.405751786e6
HOI_VSYS   = 141.5
ZOOM_VEL   = 600.0
DS_FACTOR  = 64

COLORS = {
    ('on',  'ph1', 'pol1'): '#1565c0',
    ('on',  'ph1', 'pol2'): '#42a5f5',
    ('on',  'ph2', 'pol1'): '#b71c1c',
    ('on',  'ph2', 'pol2'): '#ef9a9a',
    ('off', 'ph1', 'pol1'): '#00838f',
    ('off', 'ph1', 'pol2'): '#4db6ac',
    ('off', 'ph2', 'pol1'): '#e65100',
    ('off', 'ph2', 'pol2'): '#ffb74d',
}


def safe_path(directory, filename):
    path = os.path.join(directory, filename)
    if not os.path.exists(path):
        return path
    base, ext = os.path.splitext(filename)
    n = 1
    while os.path.exists(os.path.join(directory, f'{base}_{n}{ext}')):
        n += 1
    return os.path.join(directory, f'{base}_{n}{ext}')


def downsample(arr, factor):
    trim = (len(arr) // factor) * factor
    return arr[:trim].reshape(-1, factor).mean(axis=1)


def load_groups(fro_path):
    groups = {k: [] for k in COLORS}
    with h5py.File(fro_path, 'r') as f:
        freq_hz = f['Spectra/frequency_axis_hz'][:]
        track   = f['Quality/flag_track'][:]
        cal     = f['Quality/flag_cal'][:]
        pol     = f['Spectra/polarization'][:].astype(str)
        n       = len(track)
        spectra = f['Spectra/raw_spectrum']
        print(f'Loading {n} spectra...')
        for i in range(n):
            if i % 500 == 0:
                print(f'  {i}/{n}')
            key = (
                'on'  if track[i] == 1 else 'off',
                'ph1' if cal[i]   == 1 else 'ph2',
                pol[i]
            )
            if key in groups:
                groups[key].append(spectra[i, :])
    means  = {}
    counts = {}
    for k, specs in groups.items():
        if specs:
            means[k]  = np.mean(np.array(specs, dtype=np.float32), axis=0)
            counts[k] = len(specs)
        else:
            means[k]  = None
            counts[k] = 0
    return freq_hz, means, counts


def make_html(title, note, traces, xaxis_title, yaxis_title, vlines=None):
    traces_j = json.dumps(traces)
    vlines_j = json.dumps(vlines or [])
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: #fff; font-family: Arial, sans-serif; color: #222; }}
  h2 {{ font-size: 13px; text-align: center; padding: 6px 0 1px; color: #111; }}
  .note {{ font-size: 10px; text-align: center; color: #666; margin-bottom: 2px; }}
  #readout {{
    position: fixed; top: 52px; left: 80px;
    background: rgba(255,255,255,0.97); border: 1px solid #aaa;
    border-radius: 4px; padding: 4px 12px;
    font-size: 12px; font-family: monospace; color: #111;
    pointer-events: none; min-width: 260px; z-index: 100;
    display: none; box-shadow: 1px 1px 4px rgba(0,0,0,0.15);
  }}
  #resetbtn {{
    position: fixed; top: 8px; right: 12px;
    background: #eee; border: 1px solid #aaa; border-radius: 3px;
    padding: 3px 10px; font-size: 11px; cursor: pointer; z-index: 100;
  }}
  #resetbtn:hover {{ background: #ddd; }}
  canvas {{ display: block; cursor: crosshair; }}
  @media print {{
    #resetbtn, #readout {{ display: none !important; }}
  }}
</style>
</head>
<body>
<h2>{title}</h2>
<div class="note">{note}</div>
<div id="readout"></div>
<button id="resetbtn" onclick="resetView()">Reset zoom</button>
<canvas id="c"></canvas>
<script>
const traces = {traces_j};
const vlines = {vlines_j};
const xLabel = '{xaxis_title}';
const yLabel = '{yaxis_title}';

const PAD = {{ l:72, r:24, t:10, b:52 }};
const canvas  = document.getElementById('c');
const ctx     = canvas.getContext('2d');
const readout = document.getElementById('readout');

let view  = null;
let mouseX = null, mouseY = null;
let selStart = null, selEnd = null, selecting = false;

function fullRange() {{
  let xmin=Infinity, xmax=-Infinity, ymin=Infinity, ymax=-Infinity;
  traces.forEach(t => {{
    t.x.forEach(v => {{ if(v<xmin) xmin=v; if(v>xmax) xmax=v; }});
    t.y.forEach(v => {{ if(v<ymin) ymin=v; if(v>ymax) ymax=v; }});
  }});
  const ym=(ymax-ymin)*0.05;
  return {{xmin, xmax, ymin:ymin-ym, ymax:ymax+ym}};
}}

function getView() {{ return view || fullRange(); }}
function resetView() {{ view=null; draw(mouseX,mouseY); }}

function toCanvas(xv, yv, r) {{
  const W=canvas.width-PAD.l-PAD.r, H=canvas.height-PAD.t-PAD.b;
  return [
    PAD.l+(xv-r.xmin)/(r.xmax-r.xmin)*W,
    PAD.t+(1-(yv-r.ymin)/(r.ymax-r.ymin))*H
  ];
}}

function fromCanvas(cx, cy, r) {{
  const W=canvas.width-PAD.l-PAD.r, H=canvas.height-PAD.t-PAD.b;
  return [
    r.xmin+(cx-PAD.l)/W*(r.xmax-r.xmin),
    r.ymin+(1-(cy-PAD.t)/H)*(r.ymax-r.ymin)
  ];
}}

function inPlotArea(mx, my) {{
  const W=canvas.width-PAD.l-PAD.r, H=canvas.height-PAD.t-PAD.b;
  return mx>=PAD.l&&mx<=PAD.l+W&&my>=PAD.t&&my<=PAD.t+H;
}}

function draw(mx, my) {{
  const r=getView();
  const W=canvas.width-PAD.l-PAD.r, H=canvas.height-PAD.t-PAD.b;

  ctx.clearRect(0,0,canvas.width,canvas.height);
  ctx.fillStyle='#fff'; ctx.fillRect(0,0,canvas.width,canvas.height);
  ctx.fillStyle='#fafafa'; ctx.fillRect(PAD.l,PAD.t,W,H);

  // grid
  const nxT=8, nyT=6;
  ctx.strokeStyle='#e0e0e0'; ctx.lineWidth=1;
  for(let i=0;i<=nxT;i++) {{
    const xv=r.xmin+i*(r.xmax-r.xmin)/nxT;
    const [cx]=toCanvas(xv,0,r);
    ctx.beginPath();ctx.moveTo(cx,PAD.t);ctx.lineTo(cx,PAD.t+H);ctx.stroke();
  }}
  for(let i=0;i<=nyT;i++) {{
    const yv=r.ymin+i*(r.ymax-r.ymin)/nyT;
    const [,cy]=toCanvas(0,yv,r);
    ctx.beginPath();ctx.moveTo(PAD.l,cy);ctx.lineTo(PAD.l+W,cy);ctx.stroke();
  }}

  // tick labels
  ctx.fillStyle='#444'; ctx.font='11px Arial'; ctx.textAlign='center';
  for(let i=0;i<=nxT;i++) {{
    const xv=r.xmin+i*(r.xmax-r.xmin)/nxT;
    const [cx]=toCanvas(xv,0,r);
    ctx.fillText(xv.toFixed(1),cx,PAD.t+H+16);
  }}
  ctx.textAlign='right';
  for(let i=0;i<=nyT;i++) {{
    const yv=r.ymin+i*(r.ymax-r.ymin)/nyT;
    const [,cy]=toCanvas(0,yv,r);
    ctx.fillText(yv.toFixed(1),PAD.l-6,cy+4);
  }}

  // axis titles
  ctx.fillStyle='#222'; ctx.font='12px Arial'; ctx.textAlign='center';
  ctx.fillText(xLabel,PAD.l+W/2,canvas.height-8);
  ctx.save();
  ctx.translate(14,PAD.t+H/2); ctx.rotate(-Math.PI/2);
  ctx.fillText(yLabel,0,0); ctx.restore();

  // vlines
  vlines.forEach(vl => {{
    const [cx]=toCanvas(vl.x,0,r);
    if(cx<PAD.l||cx>PAD.l+W) return;
    ctx.strokeStyle=vl.color; ctx.lineWidth=1.2; ctx.setLineDash([6,4]);
    ctx.beginPath();ctx.moveTo(cx,PAD.t);ctx.lineTo(cx,PAD.t+H);ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle=vl.color; ctx.font='9px Arial'; ctx.textAlign='center';
    ctx.fillText(vl.label,cx,PAD.t+H+30);
  }});

  // traces (clip to plot area)
  ctx.save();
  ctx.beginPath();
  ctx.rect(PAD.l,PAD.t,W,H);
  ctx.clip();
  traces.forEach(t => {{
    ctx.strokeStyle=t.color; ctx.lineWidth=1.0;
    ctx.setLineDash(t.dash?[4,3]:[]);
    ctx.beginPath();
    let first=true;
    for(let i=0;i<t.x.length;i++) {{
      const [cx,cy]=toCanvas(t.x[i],t.y[i],r);
      if(first){{ctx.moveTo(cx,cy);first=false;}}else ctx.lineTo(cx,cy);
    }}
    ctx.stroke(); ctx.setLineDash([]);
  }});
  ctx.restore();

  // legend
  const colW=230, perRow=2;
  const ly=PAD.t+14;
  traces.forEach((t,i) => {{
    const col=i%perRow, row=Math.floor(i/perRow);
    const lx=PAD.l+10+col*colW, lyi=ly+row*16;
    ctx.strokeStyle=t.color; ctx.lineWidth=2;
    ctx.setLineDash(t.dash?[4,3]:[]);
    ctx.beginPath();ctx.moveTo(lx,lyi);ctx.lineTo(lx+20,lyi);ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle='#333'; ctx.textAlign='left'; ctx.font='10px Arial';
    ctx.fillText(t.name,lx+24,lyi+4);
  }});

  // border
  ctx.strokeStyle='#bbb'; ctx.lineWidth=1;
  ctx.strokeRect(PAD.l,PAD.t,W,H);

  // rubber-band selection rectangle
  if(selecting && selStart && selEnd) {{
    const sx=Math.min(selStart.cx,selEnd.cx);
    const sy=Math.min(selStart.cy,selEnd.cy);
    const sw=Math.abs(selEnd.cx-selStart.cx);
    const sh=Math.abs(selEnd.cy-selStart.cy);
    ctx.fillStyle='rgba(66,133,244,0.08)';
    ctx.fillRect(sx,sy,sw,sh);
    ctx.strokeStyle='rgba(66,133,244,0.8)'; ctx.lineWidth=1.5;
    ctx.setLineDash([4,3]);
    ctx.strokeRect(sx,sy,sw,sh);
    ctx.setLineDash([]);
  }}

  // crosshair + readout (only when not selecting)
  if(!selecting && mx!==null && inPlotArea(mx,my)) {{
    const [xv]=fromCanvas(mx,my,r);
    const [scx]=toCanvas(xv,0,r);
    ctx.save();
    ctx.beginPath(); ctx.rect(PAD.l,PAD.t,W,H); ctx.clip();
    ctx.strokeStyle='#bbb'; ctx.lineWidth=1; ctx.setLineDash([4,3]);
    ctx.beginPath();ctx.moveTo(scx,PAD.t);ctx.lineTo(scx,PAD.t+H);ctx.stroke();
    traces.forEach(t=>{{
      let bi=0,bd=Infinity;
      for(let i=0;i<t.x.length;i++){{const d=Math.abs(t.x[i]-xv);if(d<bd){{bd=d;bi=i;}}}}
      const [,scy]=toCanvas(xv,t.y[bi],r);
      ctx.strokeStyle=t.color;
      ctx.beginPath();ctx.moveTo(PAD.l,scy);ctx.lineTo(PAD.l+W,scy);ctx.stroke();
    }});
    ctx.setLineDash([]); ctx.restore();
    let html=`<b>${{xv.toFixed(3)}} ${{xLabel}}</b>`;
    traces.forEach(t=>{{
      let bi=0,bd=Infinity;
      for(let i=0;i<t.x.length;i++){{const d=Math.abs(t.x[i]-xv);if(d<bd){{bd=d;bi=i;}}}}
      html+=`<br><span style="color:${{t.color}}"><b>${{t.name}}</b></span>: ${{t.y[bi].toFixed(2)}} K`;
    }});
    readout.innerHTML=html; readout.style.display='block';
  }} else if(!selecting) {{
    readout.style.display='none';
  }}
}}

// Mouse events
canvas.addEventListener('mousedown', e => {{
  const rect=canvas.getBoundingClientRect();
  const mx=e.clientX-rect.left, my=e.clientY-rect.top;
  if(inPlotArea(mx,my)) {{
    selecting=true;
    selStart={{cx:mx,cy:my}};
    selEnd={{cx:mx,cy:my}};
    canvas.style.cursor='crosshair';
  }}
}});

canvas.addEventListener('mousemove', e => {{
  const rect=canvas.getBoundingClientRect();
  mouseX=e.clientX-rect.left; mouseY=e.clientY-rect.top;
  if(selecting) {{ selEnd={{cx:mouseX,cy:mouseY}}; }}
  draw(mouseX,mouseY);
}});

canvas.addEventListener('mouseup', e => {{
  if(!selecting) return;
  selecting=false;
  const rect=canvas.getBoundingClientRect();
  const mx=e.clientX-rect.left, my=e.clientY-rect.top;
  selEnd={{cx:mx,cy:my}};
  // Apply zoom only if selection is large enough (>5px)
  const dx=Math.abs(selEnd.cx-selStart.cx);
  const dy=Math.abs(selEnd.cy-selStart.cy);
  if(dx>5 && dy>5) {{
    const r=getView();
    const [x0,y0]=fromCanvas(selStart.cx,selStart.cy,r);
    const [x1,y1]=fromCanvas(selEnd.cx,selEnd.cy,r);
    view={{
      xmin:Math.min(x0,x1), xmax:Math.max(x0,x1),
      ymin:Math.min(y0,y1), ymax:Math.max(y0,y1)
    }};
  }}
  selStart=null; selEnd=null;
  draw(mouseX,mouseY);
}});

canvas.addEventListener('mouseleave', () => {{
  selecting=false; selStart=null; selEnd=null;
  mouseX=null; mouseY=null;
  readout.style.display='none';
  draw();
}});

canvas.addEventListener('dblclick', () => resetView());

function resize() {{
  canvas.width=window.innerWidth;
  canvas.height=window.innerHeight-canvas.getBoundingClientRect().top;
  draw(mouseX,mouseY);
}}
window.addEventListener('resize', resize);
resize();
</script>
</body>
</html>"""


def build_traces_fullband(freq_hz, means, counts):
    freq_flip = freq_hz[::-1]
    freq_mhz  = downsample(freq_flip, DS_FACTOR) / 1e6
    traces = []
    for k, mean in means.items():
        if mean is None:
            continue
        spec_ds = downsample(mean[::-1], DS_FACTOR)
        on_off  = 'ON'  if k[0] == 'on'  else 'OFF'
        ph      = 'Ph1' if k[1] == 'ph1' else 'Ph2'
        p       = k[2]
        traces.append({
            'x':    freq_mhz.tolist(),
            'y':    spec_ds.tolist(),
            'name': f'{ph}-{on_off} ({p}, n={counts[k]})',
            'color': COLORS[k],
            'dash':  (p == 'pol2'),
        })
    return traces


def build_traces_zoom(freq_hz, means, counts):
    freq_flip = freq_hz[::-1]
    vel_kms   = C_KMS * (1.0 - freq_flip / HI_REST_HZ)
    mask      = (vel_kms >= HOI_VSYS - ZOOM_VEL) & (vel_kms <= HOI_VSYS + ZOOM_VEL)
    vel_zoom  = vel_kms[mask]
    n_zoom    = int(np.sum(mask))
    ds        = max(1, n_zoom // 4000)
    traces = []
    for k, mean in means.items():
        if mean is None:
            continue
        spec_zoom = mean[::-1][mask]
        spec_ds   = downsample(spec_zoom, ds) if ds > 1 else spec_zoom
        vel_ds    = downsample(vel_zoom,  ds) if ds > 1 else vel_zoom
        on_off    = 'ON'  if k[0] == 'on'  else 'OFF'
        ph        = 'Ph1' if k[1] == 'ph1' else 'Ph2'
        p         = k[2]
        traces.append({
            'x':    vel_ds.tolist(),
            'y':    spec_ds.tolist(),
            'name': f'{ph}-{on_off} ({p}, n={counts[k]})',
            'color': COLORS[k],
            'dash':  (p == 'pol2'),
        })
    return traces


def main():
    fro_path = os.path.expanduser(
        sys.argv[1] if len(sys.argv) > 1 else FRO_PATH)
    out_dir  = os.path.expanduser(
        sys.argv[2] if len(sys.argv) > 2 else PLOTS_DIR)
    os.makedirs(out_dir, exist_ok=True)

    print(f'Plotting: {os.path.basename(fro_path)}')
    freq_hz, means, counts = load_groups(fro_path)

    tag = 'EFFBG_2026-06-23_6723_02-26_P217mm-EDD_bb1'

    print('Generating Plot 1: fullband...')
    traces1 = build_traces_fullband(freq_hz, means, counts)
    html1 = make_html(
        'Effelsberg 100m \u2014 HOI \u2014 P217mm-EDD \u2014 Full band, PHASE/ON-OFF',
        'T\u2090 (K, Tsys\u2248300K, not calibrated) \u2014 drag to zoom, dblclick to reset',
        traces1, 'Frequency (MHz)', 'T\u2090 (K)',
        vlines=[{'x': HI_REST_HZ/1e6, 'color': '#2e7d32',
                 'label': 'HI 1420.406 MHz'}],
    )
    out1 = safe_path(out_dir, f'{tag}_fullband_PHASE_ONOFF.html')
    with open(out1, 'w') as fh: fh.write(html1)
    print(f'Written: {out1}')

    print('Generating Plot 2: HI zoom...')
    traces2 = build_traces_zoom(freq_hz, means, counts)
    html2 = make_html(
        f'Effelsberg 100m \u2014 HOI \u2014 P217mm-EDD \u2014 HI zoom (v={HOI_VSYS}\u00b1{ZOOM_VEL:.0f} km/s)',
        'T\u2090 (K, Tsys\u2248300K, not calibrated) \u2014 drag to zoom, dblclick to reset',
        traces2, 'LSR Velocity (km/s)', 'T\u2090 (K)',
        vlines=[
            {'x': 0.0,      'color': '#2e7d32', 'label': 'v=0 km/s'},
            {'x': HOI_VSYS, 'color': '#6a1b9a', 'label': f'HOI Vsys={HOI_VSYS} km/s'},
        ],
    )
    out2 = safe_path(out_dir, f'{tag}_HI_zoom_PHASE_ONOFF.html')
    with open(out2, 'w') as fh: fh.write(html2)
    print(f'Written: {out2}')

    for p in [out1, out2]:
        subprocess.Popen(['xdg-open', p],
                         stdout=subprocess.DEVNULL,
                         stderr=subprocess.DEVNULL)
    print('\nDone.')


if __name__ == '__main__':
    main()
