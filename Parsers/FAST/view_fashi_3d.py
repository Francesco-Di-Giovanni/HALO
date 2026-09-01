#!/usr/bin/python3
# Copyright (c) 2026 Francesco Di Giovanni — HALO Project
# Released under MIT License
# https://github.com/Francesco-Di-Giovanni/HALO
#
# view_fashi_3d.py — Interactive 3D viewer for the FASHI DR2 extragalactic HI catalog
# Converts (RA, Dec, distance) to Cartesian (X, Y, Z) in Mpc and renders with Plotly.
# Output: standalone HTML file, navigable in any browser.

import os
import sys
import subprocess
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
SURVEYS_DIR  = SCRIPT_DIR
CATALOG_FILE = os.path.join(SURVEYS_DIR, "Table2_FASHI_DR2_Extragalactic_Hi_Source_Catalog.csv")
PLOTS_DIR    = os.path.join(SURVEYS_DIR, "Plots_FASHI_DR2_3D")
PLOTLY_JS    = os.path.expanduser(
    "~/.local/lib/python3.14/site-packages/plotly/package_data/plotly.min.js"
)

def safe_path(path):
    """Return path with ~ expanded and normalized."""
    return os.path.normpath(os.path.expanduser(path))

# ---------------------------------------------------------------------------
# Load catalog
# ---------------------------------------------------------------------------
def load_catalog(path):
    print(f"Loading catalog: {path}")
    df = pd.read_csv(path)
    print(f"  Total sources: {len(df):,}")

    # Drop rows with missing essential columns
    before = len(df)
    df = df.dropna(subset=["ra", "dec", "distance", "mass"])
    df = df[df["distance"] > 0]
    print(f"  Valid sources (distance > 0, no NaN): {len(df):,} ({before - len(df):,} dropped)")
    return df

# ---------------------------------------------------------------------------
# Spherical -> Cartesian (Mpc)
# Sun at origin, X toward (RA=0, Dec=0), Z toward North Pole
# ---------------------------------------------------------------------------
def to_cartesian(ra_deg, dec_deg, dist_mpc):
    ra  = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    x = dist_mpc * np.cos(dec) * np.cos(ra)
    y = dist_mpc * np.cos(dec) * np.sin(ra)
    z = dist_mpc * np.sin(dec)
    return x, y, z

# ---------------------------------------------------------------------------
# Build Plotly figure
# ---------------------------------------------------------------------------
def build_figure(df):
    x, y, z = to_cartesian(df["ra"].values, df["dec"].values, df["distance"].values)

    mass     = df["mass"].values          # log10(M_HI / M_sun)
    snr      = df["snr"].values
    name     = df["name"].values
    dist     = df["distance"].values
    v_opt    = df["v_opt"].values

    # Point size: small and uniform (156k points — keep it light)
    marker_size = 1.5

    # Tooltip
    hover = [
        f"<b>{n}</b><br>"
        f"RA={ra:.4f}°  Dec={dec:.4f}°<br>"
        f"Distance={d:.1f} Mpc<br>"
        f"log(M_HI/M☉)={m:.2f}<br>"
        f"v_opt={v:.0f} km/s<br>"
        f"SNR={s:.1f}"
        for n, ra, dec, d, m, v, s in zip(
            name, df["ra"].values, df["dec"].values,
            dist, mass, v_opt, snr
        )
    ]

    scatter = go.Scatter3d(
        x=x, y=y, z=z,
        mode="markers",
        marker=dict(
            size=marker_size,
            color=mass,
            colorscale="Viridis",
            cmin=7.0,
            cmax=11.0,
            colorbar=dict(
                title=dict(text="log(M<sub>HI</sub>/M<sub>☉</sub>)", side="right"),
                thickness=15,
                len=0.6,
            ),
            opacity=0.7,
        ),
        text=hover,
        hovertemplate="%{text}<extra></extra>",
        name="FASHI DR2",
    )

    # Observer marker at origin
    observer = go.Scatter3d(
        x=[0], y=[0], z=[0],
        mode="markers+text",
        marker=dict(size=5, color="red", symbol="diamond"),
        text=["Milky Way"],
        textposition="top center",
        textfont=dict(color="white", size=10),
        hovertemplate="Observer (Milky Way)<extra></extra>",
        name="Observer",
    )

    fig = go.Figure(data=[scatter, observer])

    fig.update_layout(
        title=dict(
            text=f"FASHI DR2 — {len(df):,} extragalactic HI sources",
            font=dict(color="white", size=16),
            x=0.5,
        ),
        paper_bgcolor="black",
        scene=dict(
            bgcolor="black",
            xaxis=dict(backgroundcolor="black", color="white", gridcolor="#333333", zerolinecolor="#555555", title="X (Mpc)"),
            yaxis=dict(backgroundcolor="black", color="white", gridcolor="#333333", zerolinecolor="#555555", title="Y (Mpc)"),
            zaxis=dict(backgroundcolor="black", color="white", gridcolor="#333333", zerolinecolor="#555555", title="Z (Mpc)"),
            
            aspectmode="cube",
        ),
        legend=dict(
            font=dict(color="white"),
            bgcolor="rgba(0,0,0,0.5)",
        ),
        margin=dict(l=0, r=0, t=40, b=0),
    )

    return fig

# ---------------------------------------------------------------------------
# Write standalone HTML with local Plotly.js
# ---------------------------------------------------------------------------
def write_html(fig, out_path, plotly_js_path):
    fig_html = fig.to_html(
        include_plotlyjs=False,
        full_html=False,
        div_id="fashi3d",
    )

    with open(plotly_js_path, "r", encoding="utf-8") as f:
        plotly_js_content = f.read()

    html = f"""<!DOCTYPE html>
<html>
<head>
<meta charset="utf-8">
<title>FASHI DR2 — 3D HI Universe</title>
<style>
  * {{ margin: 0; padding: 0; }} html, body {{ width: 100%; height: 100%; background: black; overflow: hidden; }} .plotly-graph-div {{ width: 100vw !important; height: 100vh !important; }}
  #fashi3d {{ width: 100vw; height: 100vh; }}
</style>
</head>
<body>
<script>{plotly_js_content}</script>
{fig_html}
</body>
</html>"""

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html)

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    catalog_path = safe_path(CATALOG_FILE)
    if not os.path.isfile(catalog_path):
        print(f"ERROR: catalog not found: {catalog_path}")
        sys.exit(1)

    plotly_js = safe_path(PLOTLY_JS)
    if not os.path.isfile(plotly_js):
        print(f"ERROR: plotly.min.js not found: {plotly_js}")
        sys.exit(1)

    os.makedirs(safe_path(PLOTS_DIR), exist_ok=True)

    df  = load_catalog(catalog_path)
    fig = build_figure(df)

    out_html = safe_path(os.path.join(PLOTS_DIR, "FASHI_DR2_3D_universe.html"))
    print(f"Writing HTML: {out_html}")
    write_html(fig, out_html, plotly_js)
    print("Done.")

    subprocess.Popen(["xdg-open", out_html])

if __name__ == "__main__":
    main()
