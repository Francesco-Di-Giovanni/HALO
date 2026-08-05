# view_salsa_3d.py
# 3D interactive viewer for SALSA survey data
# FRO/HALO Project — Francesco Di Giovanni, Bolzano

import numpy as np
import h5py
import pyvista as pv
from scipy.ndimage import gaussian_filter1d
import os

FRO = os.path.expanduser('~/FRO/FRO_System/SALSA_20260707_session.fro')

with h5py.File(FRO, 'r') as f:
    spectra = f['Spectra/raw_spectrum'][:]
    freq    = f['Spectra/frequency_axis_hz'][:]
    glon    = f['Pointing/glon_deg'][:]

f0 = 1420405750.0
c  = 299792.458
vlsr = (f0 - freq) / f0 * c

idx = np.argsort(glon)
glon = glon[idx]
spectra = spectra[idx]

spectra_s = gaussian_filter1d(spectra.astype(float), sigma=5, axis=1)
spectra_s -= spectra_s.min(axis=1, keepdims=True)

points = []
values = []
for i in range(len(glon)):
    for j in range(len(vlsr)):
        tb = float(spectra_s[i, j])
        if tb > 0.3:
            points.append([float(glon[i]), float(vlsr[j]), tb])
            values.append(tb)

points = np.array(points)
values = np.array(values)

cloud = pv.PolyData(points)
cloud['T_B'] = values

pl = pv.Plotter()
pl.set_background('white')
pl.add_points(cloud, scalars='T_B', cmap='jet',
              point_size=3, render_points_as_spheres=True)
pl.add_axes(xlabel='GLon', ylabel='v_LSR', zlabel='T_B')
pl.show(title='SALSA survey - 3D spectra - HALO FRO')
