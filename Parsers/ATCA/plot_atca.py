# plot_atca.py
# Interactive HTML viewer for NGC 253 HI spectrum — ATCA C1025
# Reads calibrated spectrum text file and generates a standalone Plotly.js HTML.
# Output: ~/FRO/Parsers/ATCA/Surveys/NGC253/Plots/ngc253_spectrum_interactive.html
# FRO/HALO Project — Francesco Di Giovanni (Claude AI assisted), Bolzano, Italy

import os
import sys
import json
import time
import subprocess
import numpy as np

# ── Paths ─────────────────────────────────────────────────────────────────────
_base       = os.path.expanduser('~/FRO/Parsers/ATCA/')
SPEC_TXT    = os.path.join(_base, 'Surveys/NGC253/ngc253_spectrum.txt')
PLOTLY_JS   = os.path.join(os.path.expanduser('~/FRO/Parsers/APEX/'), 'plotly-2.27.0.min.js')
OUTPUT_HTML = os.path.join(_base, 'Surveys/NGC253/Plots/ngc253_spectrum_interactive.html')
SURVEYS_DIR = os.path.join(_base, 'Surveys/NGC253/Plots/')
# ─────────────────────────────────────────────────────────────────────────────

C_KMS     = 299792.458
F_HI_MHZ  = 1420.405752
V_SYS     = 243.0


def make_html(spec_txt, plotly_js_path, output_path):
    data      = np.loadtxt(spec_txt)
    freq_mhz  = data[:, 0]
    amplitude = data[:, 1]
    vel = (C_KMS * (F_HI_MHZ - freq_mhz) / F_HI_MHZ).tolist()
    amp = amplitude.tolist()

    with open(plotly_js_path, 'r') as f:
        plotly_js = f.read()

    html = """<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>NGC 253 HI spectrum — ATCA C1025</title>
<script>""" + plotly_js + """</script>
<style>
body { margin: 0; padding: 0; font-family: Arial, sans-serif; background: #fff; display: flex; height: 100vh; overflow: hidden; }
#infobox { width: 200px; min-width: 200px; padding: 12px; border-right: 1px solid #ccc; font-size: 12px; box-sizing: border-box; overflow: hidden; }
#infobox h3 { margin: 0 0 10px 0; font-size: 13px; color: #333; }
.val-row { margin: 6px 0; }
.val-label { color: #888; font-size: 10px; text-transform: uppercase; }
.val-num { font-size: 15px; font-weight: bold; color: #1a1aff; }
hr { border: none; border-top: 1px solid #ddd; margin: 10px 0; }
#plot { flex: 1; height: 100%; min-width: 0; }
</style>
</head>
<body>
<div id="infobox">
  <h3>Cursor values</h3>
  <div class="val-row"><div class="val-label">LSRK velocity</div><div class="val-num" id="v-vel">—</div></div>
  <div class="val-row"><div class="val-label">Amplitude</div><div class="val-num" id="v-amp">—</div></div>
  <div class="val-row"><div class="val-label">Channel</div><div class="val-num" id="v-chan">—</div></div>
  <hr>
  <div class="val-row"><div class="val-label">Dataset</div><div style="font-size:10px;">ATCA C1025<br>NGC 253 HI<br>2002-02-09</div></div>
  <div class="val-row"><div class="val-label">SpW</div><div style="font-size:10px;">IF0 · 1415–1423 MHz<br>513 channels</div></div>
  <div class="val-row"><div class="val-label">v_sys</div><div style="font-size:10px;">243 km/s (LSRK)</div></div>
</div>
<div id="plot"></div>
<script>
var vel = """ + json.dumps(vel) + """;
var spec = """ + json.dumps(amp) + """;
var trace = {
    x: vel, y: spec, type: 'scatter', mode: 'lines',
    line: {color: 'blue', width: 1},
    name: 'HI amplitude', hoverinfo: 'none'
};
var layout = {
    title: {text: 'NGC 253 HI spectrum — ATCA C1025 (calibrated)', font: {size: 15}},
    xaxis: {title: 'LSRK velocity (km/s)', showspikes: true, spikemode: 'across', spikesnap: 'cursor', spikecolor: '#444', spikethickness: 1, spikedash: 'solid'},
    yaxis: {title: 'Amplitude (Jy)', showspikes: true, spikemode: 'across', spikecolor: '#ccc', spikethickness: 1, spikedash: 'dot'},
    shapes: [{type: 'line', x0: """ + str(V_SYS) + """, x1: """ + str(V_SYS) + """, y0: 0, y1: 1, yref: 'paper', line: {color: 'red', dash: 'dash', width: 1}}],
    annotations: [{x: """ + str(V_SYS) + """, y: 0.98, yref: 'paper', text: 'v_sys=""" + str(int(V_SYS)) + """ km/s', showarrow: false, xanchor: 'left', font: {color: 'red', size: 11}}],
    margin: {l: 60, r: 20, t: 50, b: 50},
    hovermode: 'x',
    dragmode: 'zoom'
};
Plotly.newPlot('plot', [trace], layout, {responsive: true, scrollZoom: false});
document.getElementById('plot').on('plotly_hover', function(data) {
    var pt = data.points[0];
    document.getElementById('v-vel').textContent = pt.x.toFixed(1) + ' km/s';
    document.getElementById('v-amp').textContent = pt.y.toExponential(3) + ' Jy';
    document.getElementById('v-chan').textContent = pt.pointIndex;
});
document.getElementById('plot').on('plotly_unhover', function() {
    document.getElementById('v-vel').textContent = '—';
    document.getElementById('v-amp').textContent = '—';
    document.getElementById('v-chan').textContent = '—';
});
</script>
</body>
</html>"""

    with open(output_path, 'w') as f:
        f.write(html)
    print('Saved:', output_path)


if __name__ == '__main__':
    if not os.path.exists(SPEC_TXT):
        print('ERROR: spectrum file not found:', SPEC_TXT)
        sys.exit(1)
    if not os.path.exists(PLOTLY_JS):
        print('ERROR: plotly-2.27.0.min.js not found:', PLOTLY_JS)
        sys.exit(1)

    make_html(SPEC_TXT, PLOTLY_JS, OUTPUT_HTML)

    # Start local HTTP server and open browser
    server = subprocess.Popen(
        ['python3', '-m', 'http.server', '8080'],
        cwd=SURVEYS_DIR,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL
    )
    time.sleep(1)
    subprocess.Popen(['xdg-open', 'http://localhost:8080/ngc253_spectrum_interactive.html'])
    print('HTTP server started on port 8080. Press Ctrl+C to stop.')
    try:
        server.wait()
    except KeyboardInterrupt:
        server.terminate()
        print('Server stopped.')
