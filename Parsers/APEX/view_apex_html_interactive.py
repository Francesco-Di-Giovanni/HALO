#!/usr/bin/env python3
"""
view_apex_html_interactive.py — APEX MBFITS 12CO(3-2) HTML interactive spectrum generator
FRO/HALO Project — Francesco Di Giovanni (Claude AI assisted)

Generates standalone HTML files using HTML5 canvas with crosshair cursor.
No external libraries required except plotly-2.27.0.min.js (already local).
Usage: python3 view_apex_html_interactive.py
"""

import os, h5py, numpy as np, json, subprocess

REST_FREQ_HZ = 345.7959899e9
C_MPS        = 299792458.0

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
FRO_DIR    = os.path.join(SCRIPT_DIR, '..', 'Surveys',
             'MBFITS_Obs_12CO_345GHz', 'MBFITS_FRO_12CO_345GHz')

fro_files = sorted([os.path.join(FRO_DIR, f)
                    for f in os.listdir(FRO_DIR) if f.endswith('.fro')])
bb1_path = [f for f in fro_files if 'bb1' in f][0]
bb2_path = [f for f in fro_files if 'bb2' in f][0]

def load_fro(path):
    with h5py.File(path, 'r') as f:
        return f['Spectra/frequency_axis_hz'][:], f['Spectra/raw_spectrum'][:]

freq1, spec1 = load_fro(bb1_path)
freq2, spec2 = load_fro(bb2_path)
mean1 = np.nanmean(spec1, axis=0)
mean2 = np.nanmean(spec2, axis=0)
vel2  = (REST_FREQ_HZ - freq2) / REST_FREQ_HZ * C_MPS / 1e3
mask  = np.abs(vel2) < 200.0
vel_z = vel2[mask]
sp_z  = mean2[mask]

step = 8

def make_html(title, note, traces, xaxis_title, yaxis_title, vline_x=None):
    traces_j = json.dumps(traces)
    vline_j  = json.dumps(vline_x)
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>{title}</title>
<style>
  * {{ box-sizing: border-box; margin: 0; padding: 0; }}
  body {{ background: #f4f4f4; font-family: Arial, sans-serif; }}
  h2 {{ font-size: 13px; text-align: center; padding: 6px 0 1px; }}
  .note {{ font-size: 10px; text-align: center; color: #777; margin-bottom: 4px; }}
  #readout {{
    position: fixed; top: 52px; left: 90px;
    background: rgba(255,255,255,0.93); border: 1px solid #aaa;
    border-radius: 4px; padding: 3px 10px;
    font-size: 12px; font-family: monospace;
    pointer-events: none; min-width: 200px; z-index: 100;
    display: none;
  }}
  canvas {{ display: block; cursor: crosshair; }}
</style>
</head>
<body>
<h2>{title}</h2>
<div class="note">{note}</div>
<div id="readout"></div>
<canvas id="c"></canvas>
<script>
const traces = {traces_j};
const vline  = {vline_j};
const xLabel = '{xaxis_title}';
const yLabel = '{yaxis_title}';

const PAD = {{ l:70, r:20, t:10, b:50 }};
const canvas = document.getElementById('c');
const ctx    = canvas.getContext('2d');
const readout = document.getElementById('readout');

function resize() {{
  canvas.width  = window.innerWidth;
  canvas.height = window.innerHeight - canvas.getBoundingClientRect().top;
  draw();
}}

function dataRange() {{
  let xmin=Infinity, xmax=-Infinity, ymin=0, ymax=-Infinity;
  traces.forEach(t => {{
    t.x.forEach(v => {{ if(v<xmin) xmin=v; if(v>xmax) xmax=v; }});
    t.y.forEach(v => {{ if(v>ymax) ymax=v; }});
  }});
  ymax *= 1.08;
  return {{xmin, xmax, ymin, ymax}};
}}

function toCanvas(xv, yv, r) {{
  const W = canvas.width - PAD.l - PAD.r;
  const H = canvas.height - PAD.t - PAD.b;
  const cx = PAD.l + (xv - r.xmin) / (r.xmax - r.xmin) * W;
  const cy = PAD.t + (1 - (yv - r.ymin) / (r.ymax - r.ymin)) * H;
  return [cx, cy];
}}

function fromCanvas(cx, cy, r) {{
  const W = canvas.width - PAD.l - PAD.r;
  const H = canvas.height - PAD.t - PAD.b;
  const xv = r.xmin + (cx - PAD.l) / W * (r.xmax - r.xmin);
  const yv = r.ymin + (1 - (cy - PAD.t) / H) * (r.ymax - r.ymin);
  return [xv, yv];
}}

let mouseX = null, mouseY = null;
let range = null;

function draw(mx, my) {{
  range = dataRange();
  const r = range;
  const W = canvas.width - PAD.l - PAD.r;
  const H = canvas.height - PAD.t - PAD.b;

  ctx.clearRect(0, 0, canvas.width, canvas.height);

  // background
  ctx.fillStyle = 'white';
  ctx.fillRect(PAD.l, PAD.t, W, H);

  // grid
  ctx.strokeStyle = '#eee'; ctx.lineWidth = 1;
  const nxTicks = 8, nyTicks = 6;
  for(let i=0; i<=nxTicks; i++) {{
    const xv = r.xmin + i*(r.xmax-r.xmin)/nxTicks;
    const [cx] = toCanvas(xv, 0, r);
    ctx.beginPath(); ctx.moveTo(cx, PAD.t); ctx.lineTo(cx, PAD.t+H); ctx.stroke();
  }}
  for(let i=0; i<=nyTicks; i++) {{
    const yv = r.ymin + i*(r.ymax-r.ymin)/nyTicks;
    const [,cy] = toCanvas(0, yv, r);
    ctx.beginPath(); ctx.moveTo(PAD.l, cy); ctx.lineTo(PAD.l+W, cy); ctx.stroke();
  }}

  // axes ticks and labels
  ctx.fillStyle = '#333'; ctx.font = '11px Arial'; ctx.textAlign = 'center';
  for(let i=0; i<=nxTicks; i++) {{
    const xv = r.xmin + i*(r.xmax-r.xmin)/nxTicks;
    const [cx] = toCanvas(xv, 0, r);
    ctx.fillText(xv.toFixed(2), cx, PAD.t+H+16);
  }}
  ctx.textAlign = 'right';
  for(let i=0; i<=nyTicks; i++) {{
    const yv = r.ymin + i*(r.ymax-r.ymin)/nyTicks;
    const [,cy] = toCanvas(0, yv, r);
    ctx.fillText(yv.toFixed(2), PAD.l-6, cy+4);
  }}

  // axis titles
  ctx.textAlign = 'center';
  ctx.fillText(xLabel, PAD.l + W/2, canvas.height - 8);
  ctx.save();
  ctx.translate(14, PAD.t + H/2);
  ctx.rotate(-Math.PI/2);
  ctx.fillText(yLabel, 0, 0);
  ctx.restore();

  // vline
  if(vline !== null) {{
    const [cvx] = toCanvas(vline, 0, r);
    ctx.strokeStyle = 'red'; ctx.lineWidth = 1.5;
    ctx.setLineDash([6,4]);
    ctx.beginPath(); ctx.moveTo(cvx, PAD.t); ctx.lineTo(cvx, PAD.t+H); ctx.stroke();
    ctx.setLineDash([]);
  }}

  // traces
  const colors = traces.map(t => t.color);
  traces.forEach(t => {{
    ctx.strokeStyle = t.color; ctx.lineWidth = 1.2;
    ctx.beginPath();
    let first = true;
    for(let i=0; i<t.x.length; i++) {{
      const [cx,cy] = toCanvas(t.x[i], t.y[i], r);
      if(first) {{ ctx.moveTo(cx,cy); first=false; }} else ctx.lineTo(cx,cy);
    }}
    ctx.stroke();

    // fill if requested
    if(t.fill) {{
      const [,cy0] = toCanvas(0, 0, r);
      ctx.lineTo(toCanvas(t.x[t.x.length-1], 0, r)[0], cy0);
      ctx.lineTo(toCanvas(t.x[0], 0, r)[0], cy0);
      ctx.closePath();
      ctx.fillStyle = t.fillcolor;
      ctx.fill();
    }}
  }});

  // legend
  let ly = PAD.t + 12;
  traces.forEach(t => {{
    ctx.strokeStyle = t.color; ctx.lineWidth = 2;
    ctx.beginPath(); ctx.moveTo(PAD.l+10, ly); ctx.lineTo(PAD.l+30, ly); ctx.stroke();
    ctx.fillStyle = '#333'; ctx.textAlign = 'left'; ctx.font = '11px Arial';
    ctx.fillText(t.name, PAD.l+34, ly+4);
    ly += 18;
  }});

  // border
  ctx.strokeStyle = '#ccc'; ctx.lineWidth = 1;
  ctx.strokeRect(PAD.l, PAD.t, W, H);

  // crosshair
  if(mx !== undefined && mx !== null) {{
    const inPlot = mx >= PAD.l && mx <= PAD.l+W && my >= PAD.t && my <= PAD.t+H;
    if(inPlot) {{
      // find x from mouse position
      const [xv] = fromCanvas(mx, my, r);
      const snapX = xv;
      const [scx] = toCanvas(snapX, 0, r);

      // vertical crosshair
      ctx.strokeStyle = '#444'; ctx.lineWidth = 1;
      ctx.setLineDash([4,3]);
      ctx.beginPath(); ctx.moveTo(scx, PAD.t); ctx.lineTo(scx, PAD.t+H); ctx.stroke();

      // horizontal crosshair per trace
      traces.forEach(t => {{
        let bi=0, bd=Infinity;
        t.x.forEach((v,i) => {{ const d=Math.abs(v-snapX); if(d<bd){{bd=d;bi=i;}} }});
        const snapY = t.y[bi];
        const [,scy] = toCanvas(snapX, snapY, r);
        ctx.strokeStyle = t.color;
        ctx.beginPath(); ctx.moveTo(PAD.l, scy); ctx.lineTo(PAD.l+W, scy); ctx.stroke();
      }});
      ctx.setLineDash([]);

      // readout
      let html = '<b>' + snapX.toFixed(5) + ' ' + xLabel + '</b>';
      traces.forEach(t => {{
        let bi=0, bd=Infinity;
        t.x.forEach((v,i) => {{ const d=Math.abs(v-snapX); if(d<bd){{bd=d;bi=i;}} }});
        html += '<br><span style="color:'+t.color+'"><b>'+t.name+'</b></span>: '+t.y[bi].toFixed(6)+' x10^14 K';
      }});
      readout.innerHTML = html;
      readout.style.display = 'block';
    }} else {{
      readout.style.display = 'none';
    }}
  }}
}}

canvas.addEventListener('mousemove', e => {{
  const rect = canvas.getBoundingClientRect();
  mouseX = e.clientX - rect.left;
  mouseY = e.clientY - rect.top;
  draw(mouseX, mouseY);
}});
canvas.addEventListener('mouseleave', () => {{
  mouseX = null; mouseY = null;
  readout.style.display = 'none';
  draw();
}});

window.addEventListener('resize', resize);
resize();
</script>
</body>
</html>"""

# traces for plot 1 — full band
traces1 = [
    {"x": (freq1[::step]/1e9).tolist(), 'y': (mean1[::step]/1e14).tolist(),
      'name': 'bb1', 'color': 'steelblue', 'fill': False},
    {'x': (freq2[::step]/1e9).tolist(), 'y': (mean2[::step]/1e14).tolist(),
      'name': 'bb2', 'color': 'darkorange', 'fill': False},
]

# traces for plot 2 — zoom velocity
traces2 = [
    {'x': vel_z[::step].tolist(), 'y': (sp_z[::step]/1e14).tolist(),
      'name': 'bb2', 'color': 'darkorange',
      'fill': True, 'fillcolor': 'rgba(255,165,0,0.15)'},
]

print("Generating Plot 1: full band...")
html1 = make_html(
    'APEX HET345-XFFTS2 \u2014 12CO(3-2) \u2014 Mean spectrum \u2014 full band (bb1 + bb2)',
    'Antenna temperature (K, nominal Tsys ~300 K, not calibrated) \u2014 Interactive: crosshair, hover for readout',
    traces1, 'Frequency (GHz)', 'Antenna temperature (x10^14 K)',
    vline_x=REST_FREQ_HZ/1e9
)
out1 = os.path.join(SCRIPT_DIR, 'APEX_12CO_mean_spectrum_full_band_interactive.html')
with open(out1, 'w') as f: f.write(html1)
print(f'  -> {{out1}}')

print("Generating Plot 2: zoom velocity...")
html2 = make_html(
    'APEX HET345-XFFTS2 \u2014 12CO(3-2) \u2014 Mean spectrum \u2014 zoom \u00b1200 km/s',
    'Antenna temperature (K, nominal Tsys ~300 K, not calibrated) \u2014 Interactive: crosshair, hover for readout',
    traces2, 'LSRK velocity (km/s)', 'Antenna temperature (x10^14 K)',
    vline_x=0.0
)
out2 = os.path.join(SCRIPT_DIR, 'APEX_12CO_mean_spectrum_zoom_interactive.html')
with open(out2, 'w') as f: f.write(html2)
print(f'  -> {{out2}}')

subprocess.Popen(['xdg-open', out1])
subprocess.Popen(['xdg-open', out2])
print('Done.')
