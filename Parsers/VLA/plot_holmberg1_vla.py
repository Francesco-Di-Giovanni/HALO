# plot_holmberg1_vla.py
# Plots for Holmberg I (DDO 63 / UGC 5139) VLA HI data — LITTLE THINGS survey
# Replicates Ott et al. 2001 (astro-ph/0110154): Figs 3, 4, 5, 6, 7, 8, 9, 10
# Author: Francesco Di Giovanni (Claude AI assisted) - FRO/HALO Project

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.patches import Ellipse, Rectangle
from astropy.io import fits
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from scipy.optimize import curve_fit
from scipy.ndimage import gaussian_filter
import warnings
warnings.filterwarnings('ignore')

# ── Paths ──────────────────────────────────────────────────────────────────────
BASE     = os.path.expanduser('~/FRO/Parsers/VLA/Surveys/Holmberg1')
RAW_DIR  = os.path.join(BASE, 'Raw')
PLOT_DIR = os.path.join(BASE, 'Plots')
os.makedirs(PLOT_DIR, exist_ok=True)

CUBE_FILE = os.path.join(RAW_DIR, 'DDO63_NA_ICL001.FITS')
MOM0_FILE = os.path.join(RAW_DIR, 'DDO63_NA_X0_P_R.FITS')
MOM1_FILE = os.path.join(RAW_DIR, 'DDO63_NA_XMOM1.FITS')
MOM2_FILE = os.path.join(RAW_DIR, 'DDO63_NA_XMOM2.FITS')

# ── Physical parameters (Ott et al. 2001, Table 3) ────────────────────────────
DYN_CENTER_RA    = (9 + 40/60 + 31.6/3600) * 15
DYN_CENTER_DEC   = 71 + 11/60 + 45/3600
MORPH_CENTER_RA  = (9 + 40/60 + 30/3600) * 15
MORPH_CENTER_DEC = 71 + 11/60 + 1.8/3600
VSYS_HEL   = 141.5
DIST_MPC   = 3.6
KPC_PER_ARCMIN = DIST_MPC * 1000 * np.tan(np.radians(1/60))
BEAM_A = 14.7 / 3600   # deg, natural beam major axis
BEAM_B = 12.7 / 3600   # deg, natural beam minor axis
BEAM_PA = 0.0
PA_MAJOR = 50.0        # deg, kinematic major axis PA (N through E)

# ── Column density conversion ──────────────────────────────────────────────────
BEAM_MAJ_RAD = np.radians(BEAM_A)
BEAM_MIN_RAD = np.radians(BEAM_B)
BEAM_AREA_SR = (np.pi / (4 * np.log(2))) * BEAM_MAJ_RAD * BEAM_MIN_RAD
LAMBDA_M  = 3e8 / 1.420405752e9
KB        = 1.380649e-23
TB_PER_JY = LAMBDA_M**2 / (2 * KB * BEAM_AREA_SR) * 1e-26
NHI_FACTOR = 1.823e18 * TB_PER_JY / 1000.0

print('Loading FITS files...')
hdu_cube = fits.open(CUBE_FILE)
hdu_mom0 = fits.open(MOM0_FILE)
hdu_mom1 = fits.open(MOM1_FILE)
hdu_mom2 = fits.open(MOM2_FILE)

cube_data = np.squeeze(hdu_cube[0].data).astype(np.float32)
mom0_data = np.squeeze(hdu_mom0[0].data).astype(np.float32)
mom1_data = np.squeeze(hdu_mom1[0].data).astype(np.float32)
mom2_data = np.squeeze(hdu_mom2[0].data).astype(np.float32)

hdr_cube = hdu_cube[0].header
hdr_mom0 = hdu_mom0[0].header

nv     = cube_data.shape[0]
crpix3 = hdr_cube['CRPIX3']
crval3 = hdr_cube['CRVAL3']
cdelt3 = hdr_cube['CDELT3']
vel_kms = (crval3 + (np.arange(nv) - (crpix3 - 1)) * cdelt3) / 1000.0

wcs2d    = WCS(hdr_mom0, naxis=2)
nhi_map  = mom0_data * NHI_FACTOR
mom1_kms = mom1_data / 1000.0
mom2_kms = mom2_data / 1000.0

cd1 = abs(hdr_mom0.get('CDELT1', -4.167e-4))
pix_scale_arcmin = cd1 * 60
pix_scale_kpc    = pix_scale_arcmin * KPC_PER_ARCMIN
beam_pix_a = BEAM_A * 3600 / (cd1 * 3600)
beam_pix_b = BEAM_B * 3600 / (cd1 * 3600)

pix_area_sr  = np.radians(cd1)**2
pix_per_beam = BEAM_AREA_SR / pix_area_sr

def sky_to_pix(ra_deg, dec_deg):
    sky = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg)
    x, y = wcs2d.world_to_pixel(sky)
    return float(x), float(y)

dyn_x,   dyn_y   = sky_to_pix(DYN_CENTER_RA,   DYN_CENTER_DEC)
morph_x, morph_y = sky_to_pix(MORPH_CENTER_RA, MORPH_CENTER_DEC)

# ── NHI mask (suppress noise pixels) ──────────────────────────────────────────
NHI_MIN = 5e19   # cm^-2, ~3 sigma threshold
nhi_masked  = np.where(nhi_map > NHI_MIN, nhi_map, np.nan)
mom1_masked = np.where(nhi_map > NHI_MIN, mom1_kms, np.nan)
mom2_masked = np.where(nhi_map > NHI_MIN, mom2_kms, np.nan)

# ── Crop to galaxy region (matching Ott et al. 2001 Fig. 5 extent) ────────────
# RA: 9h41m15s to 9h40m00s, Dec: 71d08m45s to 71d14m45s
r0, r1 = 416, 654
c0, c1 = 358, 603

nhi_crop  = nhi_masked[r0:r1, c0:c1]
nhi_crop_full = nhi_map[r0:r1, c0:c1]  # unmasked, for display in Fig 5a/5b
mom1_crop = mom1_masked[r0:r1, c0:c1]
mom2_crop = mom2_masked[r0:r1, c0:c1]
wcs_crop  = wcs2d.slice((slice(r0, r1), slice(c0, c1)))
dyn_xc,   dyn_yc   = dyn_x - c0, dyn_y - r0
morph_xc, morph_yc = morph_x - c0, morph_y - r0

# ── Smooth moment 1 to 20" (for velocity field contours, like paper) ──────────
smooth_sigma_pix = (20.0 / 3600) / cd1 / (2 * np.sqrt(2 * np.log(2)))
mom1_smooth = gaussian_filter(np.where(np.isfinite(mom1_crop), mom1_crop, 0),
                               sigma=smooth_sigma_pix)
mom1_smooth = np.where(np.isfinite(mom1_crop), mom1_smooth, np.nan)

# ── Helpers ────────────────────────────────────────────────────────────────────
def add_beam(ax, data, color='white'):
    bx = 0.06 * data.shape[1]
    by = 0.06 * data.shape[0]
    ax.add_patch(Ellipse((bx, by), width=beam_pix_b, height=beam_pix_a,
                         angle=BEAM_PA, color=color, fill=True, zorder=5))

def add_scalebar(ax, data, color='white'):
    kpc_pix = 1.0 / pix_scale_kpc
    bx = 0.07 * data.shape[1]
    by = 0.91 * data.shape[0]
    ax.annotate('', xy=(bx + kpc_pix, by), xytext=(bx, by),
                arrowprops=dict(arrowstyle='<->', color=color, lw=1.5))
    ax.text(bx + kpc_pix/2, by + 4, '1 kpc',
            ha='center', va='bottom', fontsize=7, color=color, fontweight='bold')

def make_radec_axes(ax, tick_fontsize=6):
    ax.coords[0].set_major_formatter('hh:mm:ss')
    ax.coords[1].set_major_formatter('dd:mm')
    ax.coords[0].set_ticks(spacing=2*u.arcmin)
    ax.coords[1].set_ticks(spacing=1*u.arcmin)
    ax.coords[0].set_axislabel('RA (J2000)', fontsize=9)
    ax.coords[1].set_axislabel('DEC (J2000)', fontsize=9)
    ax.coords[0].set_ticklabel(size=tick_fontsize)
    ax.coords[1].set_ticklabel(size=tick_fontsize)

print('Plotting...')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 4 — Global HI spectrum  (PNG + SVG)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 4: HI spectrum...')
flux_per_chan = np.nansum(cube_data * (nhi_map > NHI_MIN)[np.newaxis,:,:], axis=(1,2)) / pix_per_beam
mask4 = (vel_kms > 90) & (vel_kms < 200)

fig4, ax4 = plt.subplots(figsize=(8, 5))
ax4.plot(vel_kms[mask4], flux_per_chan[mask4], 'k-', lw=1.5, zorder=3)
ax4.scatter(vel_kms[mask4], flux_per_chan[mask4], color='steelblue', s=30, zorder=4)

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5*((x-mu)/sigma)**2)
try:
    good = mask4 & np.isfinite(flux_per_chan) & (flux_per_chan > 0)
    popt, _ = curve_fit(gaussian, vel_kms[good], flux_per_chan[good],
                        p0=[flux_per_chan[good].max(), VSYS_HEL, 10.0], maxfev=5000)
    v_fit = np.linspace(90, 200, 500)
    ax4.plot(v_fit, gaussian(v_fit, *popt), 'r--', lw=1.5,
             label=f'Gaussian: FWHM={2.355*abs(popt[2]):.1f} km/s, V$_c$={popt[1]:.1f} km/s')
    ax4.legend(fontsize=9)
except Exception as e:
    print(f'    Gaussian fit failed: {e}')

ax4.axvline(VSYS_HEL, color='gray', ls=':', lw=1, alpha=0.7)
ax4.set_xlabel('Heliocentric velocity [km s$^{-1}$]', fontsize=12)
ax4.set_ylabel('Flux density [Jy]', fontsize=12)
ax4.set_title('Holmberg I (DDO 63) — Global HI spectrum\nLITTLE THINGS VLA (natural weighting)', fontsize=11)
ax4.set_xlim(90, 200)
ax4.grid(True, alpha=0.3)
fig4.tight_layout()
fig4.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig4_HI_spectrum.png'), dpi=300, bbox_inches='tight')
plt.close(fig4)
print('    Saved Fig 4')

# ── Fig 4 HTML interactive ─────────────────────────────────────────────────────
v4   = vel_kms[mask4].tolist()
fl4  = flux_per_chan[mask4].tolist()
try:
    v_fit_list = np.linspace(90, 200, 500).tolist()
    fit_list   = gaussian(np.linspace(90,200,500), *popt).tolist()
    has_fit    = True
except:
    has_fit = False

html4 = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>HoI VLA — HI Spectrum</title>
<script src="plotly-2.27.0.min.js"></script>
</head><body style="background:#fff;margin:0;padding:10px">
<div id="plot" style="width:100%;height:90vh"></div>
<div id="infobox" style="position:fixed;top:10px;left:10px;background:rgba(255,255,255,0.85);border:1px solid #aaa;padding:6px 10px;font-size:13px;font-family:monospace;border-radius:4px;z-index:999;">Move mouse over plot</div>
<script>
var v = {v4};
var fl = {fl4};
var traces = [{{
  x: v, y: fl, mode: 'lines+markers',
  line: {{color:'steelblue', width:2}},
  marker: {{color:'steelblue', size:6}},
  name: 'HI flux'
}}"""
if has_fit:
    html4 += f""",{{
  x: {v_fit_list}, y: {fit_list}, mode: 'lines',
  line: {{color:'red', width:2, dash:'dash'}},
  name: 'Gaussian fit (FWHM={2.355*abs(popt[2]):.1f} km/s, Vc={popt[1]:.1f} km/s)'
}}"""
html4 += f"""];
var layout = {{
  title: 'Holmberg I (DDO 63) — Global HI spectrum<br>LITTLE THINGS VLA (natural weighting)',
  xaxis: {{title: 'Heliocentric velocity [km s⁻¹]', range:[90,200]}},
  yaxis: {{title: 'Flux density [Jy]'}},
  shapes: [{{type:'line', x0:{VSYS_HEL}, x1:{VSYS_HEL}, y0:0, y1:1,
             yref:'paper', line:{{color:'gray', dash:'dot', width:1}}}}],
  legend: {{x:0.05, y:0.95}},
  hovermode:'x unified', hoverlabel:{{bgcolor:'rgba(0,0,0,0)', bordercolor:'rgba(0,0,0,0)', font:{{color:'rgba(0,0,0,0)'}}}}, plot_bgcolor:'white', paper_bgcolor:'white'
}};
Plotly.newPlot('plot', traces, layout, {{responsive:true}});
document.getElementById('plot').on('plotly_hover', function(data) {{
  var pt = data.points[0];
  document.getElementById('infobox').innerHTML = '<b>V:</b> ' + pt.x.toFixed(2) + ' km/s<br><b>Flux:</b> ' + pt.y.toFixed(4) + ' Jy';
}});
</script>
</body></html>"""
with open(os.path.join(PLOT_DIR, 'HoI_VLA_Fig4_HI_spectrum.html'), 'w') as f:
    f.write(html4)
print('    Saved Fig 4 HTML')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 5a — Moment 0 (colour)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 5a: Moment 0...')
fig5a, ax5a = plt.subplots(figsize=(7,7), subplot_kw={'projection': wcs_crop})
im5a = ax5a.imshow(nhi_crop_full, origin='lower', cmap='inferno',
                    vmin=0, vmax=2.0e21, aspect='equal')
plt.colorbar(im5a, ax=ax5a, fraction=0.046, pad=0.04).set_label('N$_{HI}$ [cm$^{-2}$]', fontsize=9)
add_beam(ax5a, nhi_crop_full)
add_scalebar(ax5a, nhi_crop_full)
make_radec_axes(ax5a)
ax5a.set_title('Holmberg I — Integrated HI (moment 0)\nLITTLE THINGS VLA', fontsize=11)
fig5a.tight_layout()
fig5a.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig5a_mom0.png'), dpi=300, bbox_inches='tight')
plt.close(fig5a)

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 5b — Moment 0 + kinematic axes + morphological center
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 5b: Moment 0 + axes...')
fig5b, ax5b = plt.subplots(figsize=(7,7), subplot_kw={'projection': wcs_crop})
im5b = ax5b.imshow(nhi_crop_full, origin='lower', cmap='inferno',
                    vmin=0, vmax=2.0e21, aspect='equal')
plt.colorbar(im5b, ax=ax5b, fraction=0.046, pad=0.04).set_label('N$_{HI}$ [cm$^{-2}$]', fontsize=9)

PA_rad = np.radians(PA_MAJOR)
L = 180
dx_maj = -np.sin(PA_rad) * L;  dy_maj = np.cos(PA_rad) * L
dx_min = -np.sin(PA_rad+np.pi/2) * L; dy_min = np.cos(PA_rad+np.pi/2) * L
ax5b.plot([dyn_xc-dx_maj, dyn_xc+dx_maj], [dyn_yc-dy_maj, dyn_yc+dy_maj], 'w-', lw=1.5)
ax5b.plot([dyn_xc-dx_min, dyn_xc+dx_min], [dyn_yc-dy_min, dyn_yc+dy_min], 'w-', lw=1.5)
ax5b.text(dyn_xc+dx_maj*0.95, dyn_yc+dy_maj*0.95, 'Major axis',
          fontsize=7, color='white', ha='center', va='bottom')
ax5b.text(dyn_xc+dx_min*0.95, dyn_yc+dy_min*0.95, 'Minor axis',
          fontsize=7, color='white', ha='center', va='bottom')
ax5b.plot(morph_xc, morph_yc, 'w+', markersize=12, markeredgewidth=2, zorder=6)

# SE box: ~4'x4' around the small shell mentioned in paper
# Small SE shell center: alpha=9h40m35s, delta=71d10m00s (Ott et al. 2001, Sect. 4.6)
shell_x, shell_y = sky_to_pix((9 + 40/60 + 35/3600)*15, 71 + 10/60 + 0/3600)
shell_xc, shell_yc = shell_x - c0, shell_y - r0
box_pix = int(2.0 / pix_scale_arcmin)  # ~2 arcmin box  # 0.5 kpc in pixels
box_x   = shell_xc - box_pix // 2
box_y   = shell_yc - box_pix // 2
ax5b.add_patch(Rectangle((box_x, box_y), box_pix, box_pix,
                           fill=False, edgecolor='white', lw=1.5))

add_beam(ax5b, nhi_crop_full)
add_scalebar(ax5b, nhi_crop_full)
make_radec_axes(ax5b)
ax5b.set_title('Holmberg I — Integrated HI + kinematic axes\nLITTLE THINGS VLA', fontsize=11)
fig5b.tight_layout()
fig5b.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig5b_mom0_axes.png'), dpi=300, bbox_inches='tight')
plt.close(fig5b)
print('    Saved Fig 5a, 5b')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 6 — Radial HI profile  → HTML interactive
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 6: Radial profile (HTML)...')
ny, nx   = nhi_map.shape
y_i, x_i = np.mgrid[0:ny, 0:nx]
dist_as  = np.sqrt((x_i-morph_x)**2 + (y_i-morph_y)**2) * pix_scale_arcmin * 60
bin_edges   = np.arange(0, 213, 3)
bin_centers = 0.5*(bin_edges[:-1]+bin_edges[1:])
profile = []
for r0b, r1b in zip(bin_edges[:-1], bin_edges[1:]):
    m = (dist_as>=r0b)&(dist_as<r1b)&np.isfinite(nhi_masked)
    profile.append(float(np.nanmedian(nhi_masked[m])/1e21) if m.sum()>0 else float('nan'))

rad_as  = bin_centers.tolist()
rad_kpc = (bin_centers * pix_scale_arcmin / 60 * KPC_PER_ARCMIN * 60).tolist()
prof    = profile

html6 = f"""<!DOCTYPE html><html><head><meta charset="utf-8">
<title>HoI VLA — Radial HI Profile</title>
<script src="plotly-2.27.0.min.js"></script>
</head><body style="background:#fff;margin:0;padding:10px">
<div id="infobox" style="position:fixed;top:10px;left:10px;background:rgba(255,255,255,0.85);border:1px solid #aaa;padding:6px 10px;font-size:13px;font-family:monospace;border-radius:4px;z-index:999;">Move mouse over plot</div>
<div id="plot" style="width:100%;height:90vh"></div>
<script>
var rad = {rad_as};
var prof = {prof};
var traces = [{{
  x: rad, y: prof, mode: 'lines+markers',
  line: {{color:'steelblue', width:2}},
  marker: {{color:'steelblue', size:6}},
  name: 'N<sub>HI</sub>'
}}];
var layout = {{
  title: 'Holmberg I — Azimuthally-averaged HI profile<br>(origin: morphological center)',
  xaxis: {{title: 'Radius [\\'\\']', range:[0,210]}},
  yaxis: {{title: 'Column density [10²¹ cm⁻²]', rangemode:'tozero'}},
  shapes: [
    {{type:'line', x0:0, x1:210, y0:1, y1:1, line:{{color:'black',dash:'dash',width:1}}}},
    {{type:'line', x0:40, x1:40, y0:0, y1:1.5, line:{{color:'gray',dash:'dot',width:1}}}},
    {{type:'line', x0:70, x1:70, y0:0, y1:1.5, line:{{color:'gray',dash:'dot',width:1}}}}
  ],
  annotations: [
    {{x:105, y:1.02, text:'N<sub>HI</sub>=10²¹ cm⁻² (SF threshold)',
      showarrow:false, font:{{size:10}}}}
  ],
  hovermode:'x unified', hoverlabel:{{bgcolor:'rgba(0,0,0,0)', bordercolor:'rgba(0,0,0,0)', font:{{color:'rgba(0,0,0,0)'}}}}, plot_bgcolor:'white', paper_bgcolor:'white'
}};
Plotly.newPlot('plot', traces, layout, {{responsive:true}});
document.getElementById('plot').on('plotly_hover', function(data) {{
  var pt = data.points[0];
  document.getElementById('infobox').innerHTML = '<b>R:</b> ' + pt.x.toFixed(1) + '"<br><b>N<sub>HI</sub>:</b> ' + pt.y.toFixed(3) + ' ×10²¹ cm⁻²';
}});
</script>
</body></html>"""
with open(os.path.join(PLOT_DIR, 'HoI_VLA_Fig6_radial_profile.html'), 'w') as f:
    f.write(html6)
print('    Saved Fig 6 HTML')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 7 — Velocity field (smoothed moment 1 contours on moment 0)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 7: Velocity field...')
fig7, ax7 = plt.subplots(figsize=(7,7), subplot_kw={'projection': wcs_crop})
im7 = ax7.imshow(nhi_crop, origin='lower', cmap='inferno',
                  vmin=0, vmax=2.0e21, aspect='equal')
plt.colorbar(im7, ax=ax7, fraction=0.046, pad=0.04).set_label('N$_{HI}$ [cm$^{-2}$]', fontsize=9)
v_levels = np.arange(125, 162, 5)
valid7   = np.isfinite(mom1_smooth)
if valid7.sum() > 100:
    cs = ax7.contour(np.where(valid7, mom1_smooth, np.nan),
                      levels=v_levels, colors='white', linewidths=0.9, origin='lower')
    ax7.clabel(cs, fmt='%d', fontsize=6, inline=True, colors='white')
add_beam(ax7, nhi_crop)
add_scalebar(ax7, nhi_crop)
make_radec_axes(ax7)
ax7.set_title('Holmberg I — Velocity field (moment 1, smoothed 20\'\')\nLITTLE THINGS VLA', fontsize=11)
fig7.tight_layout()
fig7.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig7_velocity_field.png'), dpi=300, bbox_inches='tight')
plt.close(fig7)
print('    Saved Fig 7')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 8 — Position-velocity cuts (major and minor axis)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 8: PV cuts...')

def extract_pv_cut(cube, center_x, center_y, pa_deg, width_pix=3):
    """Extract PV diagram along axis defined by PA through center pixel."""
    ny_c, nx_c = cube.shape[1], cube.shape[2]
    pa_rad = np.radians(pa_deg)
    # Length of cut: diagonal of image
    L = int(np.sqrt(nx_c**2 + ny_c**2) / 2)
    offsets = np.arange(-L, L+1)
    pv = np.full((cube.shape[0], len(offsets)), np.nan, dtype=np.float32)
    for i, off in enumerate(offsets):
        # pixel along axis: dx=-sin(PA), dy=cos(PA) (N through E, y=N)
        px = center_x - off * np.sin(pa_rad)
        py = center_y + off * np.cos(pa_rad)
        # average over width perpendicular to axis
        vals = []
        for dw in range(-width_pix//2, width_pix//2+1):
            wx = px + dw * np.cos(pa_rad)
            wy = py + dw * np.sin(pa_rad)
            ix, iy = int(round(wx)), int(round(wy))
            if 0 <= ix < nx_c and 0 <= iy < ny_c:
                vals.append(cube[:, iy, ix])
        if vals:
            pv[:, i] = np.nanmean(vals, axis=0)
    return offsets, pv

# Use morphological center for PV cuts (as in paper Fig. 8)
offsets_maj, pv_maj = extract_pv_cut(cube_data, morph_x, morph_y, PA_MAJOR)
offsets_min, pv_min = extract_pv_cut(cube_data, morph_x, morph_y, PA_MAJOR+90)

# Noise
noise_ch = cube_data[(vel_kms<100)|(vel_kms>190)]
rms = float(np.nanstd(noise_ch))

# Trim to meaningful range (±250 arcsec)
max_off_pix = int(250 / (pix_scale_arcmin*60))
trim = max_off_pix
off_maj_as = offsets_maj * pix_scale_arcmin * 60
off_min_as = offsets_min * pix_scale_arcmin * 60
mask_maj = np.abs(off_maj_as) <= 250
mask_min = np.abs(off_min_as) <= 250

fig8, (ax8a, ax8b) = plt.subplots(2, 1, figsize=(8, 10))
for ax, pv, off_as, label in [
        (ax8a, pv_maj, off_maj_as, 'Major Axis'),
        (ax8b, pv_min, off_min_as, 'Minor Axis')]:
    off_sel = off_as[mask_maj if label=='Major Axis' else mask_min]
    pv_sel  = pv[:, mask_maj if label=='Major Axis' else mask_min]
    # Convert Jy/beam to K
    pv_K = pv_sel * TB_PER_JY
    ax.imshow(pv_K, origin='lower', aspect='auto', cmap='inferno',
               vmin=-2*rms*TB_PER_JY, vmax=15*rms*TB_PER_JY,
               extent=[off_sel[0], off_sel[-1], vel_kms[0], vel_kms[-1]])
    ax.set_ylim(95, 195)
    ax.set_xlim(-250, 250)
    ax.axhline(VSYS_HEL, color='white', ls='--', lw=0.8, alpha=0.6)
    ax.axvline(0, color='white', ls='--', lw=0.8, alpha=0.6)
    ax.set_ylabel('V$_{hel}$ [km s$^{-1}$]', fontsize=11)
    ax.set_xlabel('Angular offset [\'\']', fontsize=11)
    ax.tick_params(labelsize=8)
    ax.text(0.05, 0.07, label, transform=ax.transAxes,
            color='white', fontsize=10, fontweight='bold')

fig8.suptitle('Holmberg I — Position-velocity cuts\nLITTLE THINGS VLA', fontsize=11)
fig8.tight_layout()
fig8.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig8_pv_cuts.png'), dpi=300, bbox_inches='tight')
plt.close(fig8)
print('    Saved Fig 8')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 9 — Rotation curve V×sin(i) (simplified from moment 1)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 9: Rotation curve...')
# Derive V*sin(i) from smoothed moment 1 along major axis
# Tilt-ring approximation: measure mean velocity in annuli on approaching/receding sides
off_maj_as_full = offsets_maj * pix_scale_arcmin * 60
off_kpc = off_maj_as_full * pix_scale_kpc / pix_scale_arcmin / 60 * pix_scale_arcmin * 60

# Extract velocity along major axis from smoothed moment 1
v_along = []
r_along = []
for i, off in enumerate(offsets_maj):
    px = morph_x - off * np.sin(np.radians(PA_MAJOR))
    py = morph_y + off * np.cos(np.radians(PA_MAJOR))
    ix, iy = int(round(px - c0)), int(round(py - r0))
    if 0<=ix<mom1_smooth.shape[1] and 0<=iy<mom1_smooth.shape[0]:
        v = mom1_smooth[iy, ix]
        if np.isfinite(v):
            r_as = off * pix_scale_arcmin * 60
            v_along.append((r_as, v))

if v_along:
    r_as_arr = np.array([x[0] for x in v_along])
    v_arr    = np.array([x[1] for x in v_along])
    # V*sin(i) = |V_obs - V_sys|
    vrot_arr = np.abs(v_arr - VSYS_HEL)
    r_kpc_arr = r_as_arr * pix_scale_kpc / (pix_scale_arcmin * 60)

    # Keep only |r| > 5 arcsec and |r| < 220 arcsec
    msk9 = (np.abs(r_as_arr) > 5) & (np.abs(r_as_arr) < 220)
    r_plot  = np.abs(r_as_arr[msk9])
    v_plot  = vrot_arr[msk9]
    rk_plot = np.abs(r_kpc_arr[msk9])

    # Sort by radius
    sort_idx = np.argsort(r_plot)
    r_plot, v_plot, rk_plot = r_plot[sort_idx], v_plot[sort_idx], rk_plot[sort_idx]

    fig9, ax9 = plt.subplots(figsize=(8, 5))
    ax9.errorbar(r_plot, v_plot, fmt='ko', markersize=5, capsize=3,
                  label='V × sin(i) (major axis)')
    ax9.set_xlabel('Radius [\'\']', fontsize=12)
    ax9.set_ylabel('V × sin(i) [km s$^{-1}$]', fontsize=12)
    ax9.set_title('Holmberg I — Observed rotation curve\nLITTLE THINGS VLA', fontsize=11)
    ax9.set_xlim(0, 220)
    ax9.set_ylim(0, max(20, v_plot.max()*1.2))
    ax9.grid(True, alpha=0.3)
    ax9.legend(fontsize=9)
    ax9_top = ax9.twiny()
    ax9_top.set_xlim(0, 220*pix_scale_kpc/(pix_scale_arcmin*60)*pix_scale_arcmin*60)
    ax9_top.set_xlabel('Radius [kpc]', fontsize=11)
    fig9.tight_layout()
    fig9.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig9_rotation_curve.png'), dpi=300, bbox_inches='tight')
    plt.close(fig9)
    print('    Saved Fig 9')
else:
    print('    Fig 9: no valid velocity data along major axis')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 10 — Velocity dispersion (moment 2, masked)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 10: Velocity dispersion...')
fig10, ax10 = plt.subplots(figsize=(7,7), subplot_kw={'projection': wcs_crop})
im10 = ax10.imshow(mom2_crop, origin='lower', cmap='plasma',
                    vmin=5, vmax=20, aspect='equal')
plt.colorbar(im10, ax=ax10, fraction=0.046, pad=0.04).set_label('Velocity dispersion [km s$^{-1}$]', fontsize=9)
add_beam(ax10, mom2_crop)
add_scalebar(ax10, mom2_crop)
make_radec_axes(ax10)
ax10.set_title('Holmberg I — HI velocity dispersion (moment 2)\nLITTLE THINGS VLA', fontsize=11)
fig10.tight_layout()
fig10.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig10_velocity_dispersion.png'), dpi=300, bbox_inches='tight')
plt.close(fig10)
print('    Saved Fig 10')

# ══════════════════════════════════════════════════════════════════════════════
# FIG. 3 — Channel maps (3×N grid, viridis)
# ══════════════════════════════════════════════════════════════════════════════
print('  Fig 3: Channel maps...')
chan_mask    = (vel_kms >= 120) & (vel_kms <= 175)
chan_indices = np.where(chan_mask)[0]
if len(chan_indices) > 12:
    step = len(chan_indices) // 12
    chan_indices = chan_indices[::step][:12]

noise_ch2 = cube_data[(vel_kms<100)|(vel_kms>190)]
rms2 = float(np.nanstd(noise_ch2))
vmin_ch, vmax_ch = -2*rms2, 8*rms2

ncols = 3
nrows = (len(chan_indices)+ncols-1)//ncols
fig3, axes3 = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows),
                             subplot_kw={'projection': wcs_crop})
axes3 = np.array(axes3).reshape(nrows, ncols)

for idx, ch in enumerate(chan_indices):
    row, col = divmod(idx, ncols)
    ax = axes3[row, col]
    chan_map = cube_data[ch, r0:r1, c0:c1]
    ax.imshow(chan_map, origin='lower', cmap='inferno',
               vmin=vmin_ch, vmax=vmax_ch, aspect='equal')
    ax.text(0.05, 0.93, f'{vel_kms[ch]:.0f} km/s',
            transform=ax.transAxes, fontsize=7, color='white', va='top')
    ax.add_patch(Ellipse((0.06*chan_map.shape[1], 0.06*chan_map.shape[0]),
                          width=beam_pix_b, height=beam_pix_a,
                          angle=BEAM_PA, color='white', fill=True, zorder=5))
    if row == nrows-1:
        ax.coords[0].set_axislabel('RA (J2000)', fontsize=6)
        ax.coords[0].set_major_formatter('hh:mm:ss')
        ax.coords[0].set_ticks(spacing=2*u.arcmin)
        ax.coords[0].set_ticklabel(size=5)
    else:
        ax.coords[0].set_ticklabel_visible(False)
        ax.coords[0].set_axislabel('')
    if col == 0:
        ax.coords[1].set_axislabel('DEC (J2000)', fontsize=6)
        ax.coords[1].set_major_formatter('dd:mm')
        ax.coords[1].set_ticks(spacing=1*u.arcmin)
        ax.coords[1].set_ticklabel(size=5)
    else:
        ax.coords[1].set_ticklabel_visible(False)
        ax.coords[1].set_axislabel('')

for idx in range(len(chan_indices), nrows*ncols):
    axes3[divmod(idx, ncols)].set_visible(False)

fig3.suptitle('Holmberg I — HI channel maps (km s$^{-1}$ heliocentric)\nLITTLE THINGS VLA (natural weighting)', fontsize=11)
fig3.tight_layout()
fig3.savefig(os.path.join(PLOT_DIR, 'HoI_VLA_Fig3_channel_maps.png'), dpi=200, bbox_inches='tight')
plt.close(fig3)
print('    Saved Fig 3')

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
hdu_cube.close(); hdu_mom0.close(); hdu_mom1.close(); hdu_mom2.close()
print('\nAll plots saved to:', PLOT_DIR)
for f in sorted(os.listdir(PLOT_DIR)):
    fpath = os.path.join(PLOT_DIR, f)
    print(f'  {f} ({os.path.getsize(fpath)//1024} KB)')
