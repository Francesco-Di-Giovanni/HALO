# view_hi4pi_surface.py
# 3D terrain-style surface viewer for HI4PI FITS cube using PyVista
# Height and color both driven by peak brightness temperature (max over velocity axis)
# per sky position — avoids volume-rendering accumulation artifacts.
# FRO/HALO Project — Francesco Di Giovanni, Bolzano
import numpy as np
from astropy.io import fits
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
import os

blue_white_cmap = LinearSegmentedColormap.from_list('blue_white', ['navy', 'dodgerblue', 'white'])

CUBE_FILE = os.path.expanduser('~/FRO/Parsers/HI4PI/CAR_D05.fits')

print('Loading HI4PI cube...')
data = fits.getdata(CUBE_FILE)  # shape: (vel, glat, glon)
data = np.nan_to_num(data, nan=0.0)
print(f'Cube shape: {data.shape}')

print('Computing peak temperature map (max over velocity axis)...')
peak_T = data.max(axis=0)  # shape: (glat, glon)
n_glat, n_glon = peak_T.shape
print(f'Peak map shape: {peak_T.shape}')
print(f'Peak T range: {peak_T.min():.2f} to {peak_T.max():.2f} K')

# Build a regular grid in pixel space for glon/glat
x = np.arange(n_glon)
y = np.arange(n_glat)
X, Y = np.meshgrid(x, y)

# Height exaggeration: peak_T is in K, pixel spacing is 1 -- scale height so
# terrain relief is visually comparable to the spatial extent.
height_exaggeration = 1.5
Z = peak_T * height_exaggeration / peak_T.max() * (min(n_glon, n_glat) * 0.3)

grid = pv.StructuredGrid(X, Y, Z)
grid['Peak T_B'] = peak_T.flatten(order='F')  # StructuredGrid uses Fortran order for point data matching X,Y,Z construction

print('Launching 3D terrain viewer...')

pl = pv.Plotter()
pl.add_mesh(
    grid,
    scalars='Peak T_B',
    cmap=blue_white_cmap,
    clim=[10.0, 60.0],
    smooth_shading=True,
    show_scalar_bar=True,
    scalar_bar_args={'title': 'Peak T_B (K)', 'color': 'white'},
)
pl.add_axes()
pl.show_grid()
pl.set_background('black')
pl.camera.zoom(1.5)

HTML_PATH = os.path.expanduser('~/FRO/Parsers/HI4PI/hi4pi_peak_terrain_interactive.html')
pl.export_html(HTML_PATH)
print(f'Interactive HTML saved: {HTML_PATH}')
print('Open it in any web browser to rotate/zoom with the mouse -- no Python needed.')

# --- SpacePilot Pro (pyspacemouse) integration, with mouse-only fallback ---
import time
try:
    import pyspacemouse
    spacemouse_ok = pyspacemouse.open()
except Exception as e:
    spacemouse_ok = False
    print(f'SpacePilot Pro not available ({e}); falling back to mouse-only navigation.')

if spacemouse_ok:
    print('SpacePilot Pro connected: move to rotate, twist to zoom, press any button to reset view.')
    print('Close the window (or Ctrl+C in terminal) to exit.')
    rot_sensitivity = 0.3   # degrees per unit per frame -- kept small, loop runs ~100x/sec
    zoom_sensitivity = 0.01
    deadzone = 0.05
    current_elevation = 0.0  # tracked ourselves, clamped to avoid gimbal-lock warning spam
    ELEV_LIMIT = 85.0
    pl.show(title='HI4PI Peak Temperature Terrain - HALO FRO', auto_close=False, interactive_update=True)
    try:
        while True:
            state = spacemouse_ok.read()
            if state is not None:
                if abs(state.yaw) > deadzone:
                    pl.camera.azimuth += state.yaw * rot_sensitivity
                if abs(state.pitch) > deadzone:
                    delta = state.pitch * rot_sensitivity
                    new_elev = max(-ELEV_LIMIT, min(ELEV_LIMIT, current_elevation + delta))
                    applied = new_elev - current_elevation
                    if abs(applied) > 1e-6:
                        pl.camera.elevation += applied
                        current_elevation = new_elev
                if abs(state.roll) > deadzone:
                    pl.camera.roll += state.roll * rot_sensitivity
                if abs(state.z) > deadzone:
                    zoom_factor = 1.0 + state.z * zoom_sensitivity
                    if zoom_factor > 0:
                        pl.camera.zoom(zoom_factor)
                if any(state.buttons):
                    pl.reset_camera()
                    pl.camera.zoom(1.5)
                    current_elevation = 0.0
            pl.update()
            time.sleep(0.01)
    except Exception:
        import traceback
        traceback.print_exc()
    finally:
        try:
            pl.close()
        except Exception:
            pass
else:
    print('Use mouse to rotate, scroll to zoom.')
    print('Close window to exit.')
    SCREENSHOT_PATH = os.path.expanduser('~/FRO/Parsers/HI4PI/hi4pi_peak_terrain_example.png')
    pl.show(title='HI4PI Peak Temperature Terrain - HALO FRO', screenshot=SCREENSHOT_PATH)
    print(f'Screenshot saved: {SCREENSHOT_PATH}')
