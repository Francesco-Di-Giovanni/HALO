# view_salsa_terrain.py
# 3D terrain-style surface viewer for the SALSA galactic HI survey
# (SALSA_Survey_Galactic_5deg_120s). Shows the baseline-subtracted peak map
# (method B: DC+RFI masked, then per-spectrum baseline subtraction removes
# passband slope/ripple before taking the peak). Height and color are both
# driven by peak antenna temperature (K, nominal, not calibrated).
# X/Y axes are real Galactic longitude/latitude in degrees. The Z axis is a
# relative, arbitrarily-scaled height for readability -- it is NOT a physical
# quantity; the real value (K, nominal) is encoded only in the color/colorbar.
# The Galactic plane (b=0 deg) is marked with a highlighted ridge line that
# follows the terrain surface exactly at that latitude.
# A flat tray sits below the terrain, extending past the near GLat edge
# toward the observer, carrying the longitude gridlines (every 10 deg) so
# they read on a clean surface, detached from the terrain itself.
# FRO/HALO Project - Francesco Di Giovanni, Bolzano

import os
import numpy as np
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap

blue_white_cmap = LinearSegmentedColormap.from_list('blue_white', ['navy', 'dodgerblue', 'white'])

NPZ_FILE = os.path.expanduser(
    '~/FRO/Parsers/SALSA/Surveys/SALSA_Survey_Galactic_5deg_120s/salsa_terrain_grid.npz'
)

print('Loading SALSA survey terrain grid...')
d = np.load(NPZ_FILE)
glon_values = d['glon_values']  # real GLon degrees, 1D
glat_values = d['glat_values']  # real GLat degrees, 1D
peak_grid_b = d['peak_grid_b']  # shape (n_glat, n_glon) -- baseline-subtracted

n_glat, n_glon = peak_grid_b.shape
print(f'Grid shape: {n_glat} GLat x {n_glon} GLon')
n_gaps = int(np.sum(np.isnan(peak_grid_b)))
if n_gaps:
    print(f'NOTE: {n_gaps} grid cells have no observation (shown as flat/zero).')

# Replace NaN gaps with 0 for display only -- these are known coverage gaps
# (e.g. Aquila-Cygnus GLon 40-80 deferred to September), not zero emission.
peak_display = np.nan_to_num(peak_grid_b, nan=0.0)

# Real Galactic coordinates as the X/Y mesh -- axes show true degrees.
X, Y = np.meshgrid(glon_values, glat_values)

vmin = 0.0
vmax = float(peak_display.max())

glon_span = float(glon_values.max() - glon_values.min())
glat_span = float(glat_values.max() - glat_values.min())
BASE_HEIGHT_SCALE = min(glon_span, glat_span) * 0.3

Z = peak_display / vmax * BASE_HEIGHT_SCALE  # scale factor = 1.0 at construction
grid = pv.StructuredGrid(X, Y, Z)
grid['Peak T_B'] = peak_display.flatten(order='F')
base_z = Z.copy()

# Galactic plane reference (b = 0 deg): a ridge line that follows the terrain
# surface exactly at GLat = 0 -- always fully visible across the whole
# longitude range, since it sits directly on the terrain.
glat_zero_idx = int(np.argmin(np.abs(glat_values - 0.0)))
row_z = base_z[glat_zero_idx, :].copy()

ridge_pts = np.column_stack([glon_values, np.zeros(n_glon), row_z])
ridge = pv.PolyData()
ridge.points = ridge_pts
ridge.lines = np.hstack([[n_glon], np.arange(n_glon)])

# Flat tray: a horizontal plane just below the terrain, extending past its
# near GLat edge toward the observer. Longitude gridlines (every 10 deg) lie
# flat on top of it -- trays are horizontal, no bend.
TRAY_Z = -0.1 * BASE_HEIGHT_SCALE
LIP_WIDTH = 0.7 * glat_span

glat_edge = float(glat_values.max())      # near edge of the terrain data
tray_far_y = glat_edge + LIP_WIDTH        # outer edge of the tray, toward observer

tray_Xg, tray_Yg = np.meshgrid(glon_values, [tray_far_y, glat_edge])
tray_Zg = np.full_like(tray_Xg, TRAY_Z)
tray = pv.StructuredGrid(tray_Xg, tray_Yg, tray_Zg)

# Longitude gridlines on the tray, every 10 deg
lon_ticks = np.arange(np.ceil(glon_values.min() / 10) * 10, glon_values.max() + 0.01, 10)
tray_grid_pts, tray_grid_lines = [], []
for i, lon in enumerate(lon_ticks):
    tray_grid_pts.append([lon, tray_far_y, TRAY_Z])
    tray_grid_pts.append([lon, glat_edge, TRAY_Z])
    tray_grid_lines.append([2, 2 * i, 2 * i + 1])
tray_gridlines = pv.PolyData()
tray_gridlines.points = np.array(tray_grid_pts)
tray_gridlines.lines = np.hstack(tray_grid_lines)

# Explicit axis bounds: GLon/GLat left unpadded so tick labels land on clean
# multiples of 10/5 deg. Z is extended below the tray, and above to stay
# valid up to the top of the height-scale slider range (5x).
AXIS_BOUNDS = (
    float(glon_values.min()), float(glon_values.max()),
    float(glat_values.min()), float(glat_values.max()),
    TRAY_Z - 0.1 * BASE_HEIGHT_SCALE, 5.2 * BASE_HEIGHT_SCALE,
)

pl = pv.Plotter()

pl.add_text('SALSA survey - masked + baseline subtracted', font_size=10, color='white')
pl.add_mesh(
    grid,
    scalars='Peak T_B',
    cmap=blue_white_cmap,
    clim=[vmin, vmax],
    smooth_shading=True,
    show_scalar_bar=True,
    scalar_bar_args={
        'title': 'Peak T (K, nominal)',
        'color': 'white',
        'vertical': True,
        'position_x': 0.08,
        'position_y': 0.15,
        'width': 0.05,
        'height': 0.65,
    },
)
pl.add_mesh(ridge, color='yellow', line_width=4, render_lines_as_tubes=True, name='galactic_plane')
pl.add_mesh(tray, color='gray', opacity=0.15, name='tray')
pl.add_mesh(tray_gridlines, color='white', line_width=1, opacity=0.5, name='tray_gridlines')
pl.show_bounds(
    location='origin',
    xtitle='Galactic longitude l (deg)',
    ytitle='Galactic latitude b (deg)',
    ztitle='Height (relative, arbitrary -- see colorbar for K)',
    color='white',
    n_xlabels=13,
    n_ylabels=9,
    bounds=AXIS_BOUNDS,
)
pl.set_background('black')
pl.view_xy()
pl.camera.zoom(1.3)


def update_height_scale(value):
    # Snap to 0.1 steps so the slider settles cleanly on round values (e.g. 1.0)
    value = round(value * 10) / 10
    grid.points[:, 2] = (base_z * value).flatten(order='F')
    ridge.points[:, 2] = row_z * value
    pl.render()


height_scale_slider = pl.add_slider_widget(
    update_height_scale,
    rng=[0.1, 5.0],
    value=1.0,
    title='Height scale (x)',
    pointa=(0.75, 0.1),
    pointb=(0.95, 0.1),
    style='modern',
    color='white',
    fmt='%.1f',
)


def reset_height_scale():
    height_scale_slider.GetRepresentation().SetValue(1.0)
    update_height_scale(1.0)
    pl.render()


pl.add_key_event('r', reset_height_scale)

print('Launching 3D terrain viewer...')
print('Use the slider to change the height scale live.')
print("Press 'r' to reset height scale to 1.0.")
print('Rotate with left-drag, zoom with scroll, pan with Shift+left-drag (or middle-drag).')
print('Yellow ridge line follows the terrain at the Galactic plane (b = 0 deg).')
print('Longitude gridlines are on the tray below the terrain.')

SCREENSHOT_PATH = os.path.expanduser(
    '~/FRO/Parsers/SALSA/salsa_terrain_example.png'
)
pl.show(
    title='SALSA Survey Peak Terrain - HALO FRO',
    screenshot=SCREENSHOT_PATH,
)
print(f'Screenshot saved: {SCREENSHOT_PATH}')
