# view_hi4pi_cube.py
# 3D interactive viewer for HI4PI FITS cube using PyVista
# SpacePilot Pro support via pyspacemouse
# FRO/HALO Project — Francesco Di Giovanni, Bolzano

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import pyvista as pv
from matplotlib.colors import LinearSegmentedColormap
blue_white_cmap = LinearSegmentedColormap.from_list('blue_white', ['navy', 'dodgerblue', 'white'])
import sys
import os

CUBE_FILE = os.path.expanduser(
    '~/FRO/Parsers/HI4PI/CAR_D05.fits')

print('Loading HI4PI cube...')
h = fits.open(CUBE_FILE)
data = h[0].data  # shape: (vel, glat, glon)
hdr = h[0].header
wcs = WCS(hdr)
h.close()

print(f'Cube shape: {data.shape}')
print(f'Data range: {np.nanmin(data):.2f} to {np.nanmax(data):.2f} K')

# Replace NaN with 0
data = np.nan_to_num(data, nan=0.0)

# Threshold: show only significant emission
threshold = 0.5  # K
data_thresh = np.where(data > threshold, data, 0.0)

print(f'Creating PyVista volume...')
grid = pv.ImageData()
grid.dimensions = np.array(data_thresh.shape[::-1]) + 1
n_vel, n_glat, n_glon = data_thresh.shape
depth_exaggeration = 8.0  # visual only: increase to make velocity-axis structure more prominent
spacing_z = (n_glon / n_vel) * depth_exaggeration
grid.spacing = (1, 1, spacing_z)
print(f'Velocity axis visual spacing: {spacing_z:.4f} (cube depth roughly matches spatial extent)')
grid.cell_data['T_B'] = data_thresh.flatten(order='C')

print('Launching 3D viewer...')
print('Use mouse to rotate, scroll to zoom.')
print('Close window to exit.')

pl = pv.Plotter()
pl.add_volume(grid, scalars='T_B',
              cmap=blue_white_cmap,
              opacity='sigmoid',
              clim=[10.0, 19.0],
              shade=False,  # disable gradient-based lighting: color should reflect data value only, not steepness
              show_scalar_bar=True,
              scalar_bar_args={'title': 'T_B (K)', 'color': 'white'})
pl.add_axes()
pl.show_grid()
pl.set_background('black')
pl.camera.zoom(1.8)  # start closer, avoids manual zoom softening details
pl.show(title='HI4PI Cube - HALO FRO')
