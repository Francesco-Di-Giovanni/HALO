#!/usr/bin/env python3
# view_salsa_terrain.py
# 3D terrain-style surface viewer for the SALSA galactic HI survey
# (e.g. SALSA_Survey_Galactic_5deg_120s). Shows the baseline-subtracted
# peak map (method B: DC+RFI masked, then per-spectrum baseline subtraction
# removes passband slope/ripple before taking the peak), read from the
# .npz grid written by build_salsa_terrain_grid.py, or built on the fly
# from a HALO .fro file with the same method B.
#
# Height and colour are both driven by the peak antenna temperature
# (K, nominal, not calibrated). X/Y are real Galactic longitude and
# latitude in degrees. The height is relative and arbitrary (x1.0: the
# highest peak rises to 30 % of the smaller angular span); the real value
# is the colour.
#
# Display options (all optional, Qt side panel):
#  - smooth surface: a cubic spline through the data nodes (exact at the
#    nodes), evaluated on a mesh SMOOTH_FACTOR times finer, for display
#    only; clipped to the data range. Unticked: the plain grid of the data
#    nodes, with linear contours;
#  - Galactic plane (b = 0, yellow ridge), graticule every 10 deg on the
#    surface, contours every CONTOUR_STEP_K kelvin, graticule on the floor,
#    axis labels just outside the four floor edges.
# Surface, ridge, graticule and contours are built on the same triangle
# mesh, so the lines lie exactly on the rendered surface.
# Dark background on purpose: the colour map ends in white (peaks).
# Output: PNG images (x1/x2/x3) and a self-contained interactive HTML page
# (Plotly, offline) that opens with the current view and options, both
# through a "Save as" window; "Preview in browser" writes one temporary
# file outside the project. In the HTML page the hover values come from
# the data nodes only. (SVG export was removed: VTK's vector exporter
# garbles the 2D texts.)
# The window reopens with the size and position of the last session.
#
# Usage: view_salsa_terrain.py <file.fro | salsa_terrain_grid.npz>
#        (normally started by the HALO Launcher)
#
# FRO/HALO Project - Francesco Di Giovanni, Bolzano

import os
import sys
import tempfile
from datetime import datetime

import numpy as np
import pyvista as pv
from scipy.interpolate import RectBivariateSpline
from matplotlib.colors import LinearSegmentedColormap
from PyQt5.QtCore import Qt, QSettings, QUrl
from PyQt5.QtGui import QDesktopServices
from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QFrame,
                             QHBoxLayout, QVBoxLayout, QLabel, QSlider,
                             QPushButton, QCheckBox, QMessageBox,
                             QComboBox, QFileDialog, QScrollArea)
from pyvistaqt import QtInteractor
from vtkmodules.vtkCommonCore import vtkDoubleArray
from vtkmodules.vtkRenderingCore import vtkTextActor

BLUE_WHITE = LinearSegmentedColormap.from_list(
    'blue_white', ['navy', 'dodgerblue', 'white'])
PLOTLY_COLORSCALE = [[0.0, 'navy'], [0.5, 'dodgerblue'], [1.0, 'white']]
FLAT_LIGHT = dict(ambient=1.0, diffuse=0.0, specular=0.0,
                  roughness=1.0, fresnel=0.0)

HEIGHT_MIN, HEIGHT_MAX, HEIGHT_STEP = 0.1, 5.0, 0.1
LABEL_FONT_SIZE = 14         # axis labels, pixels on screen
TITLE_FONT_SIZE = 20         # title, upper left
CBAR_TITLE_FONT_SIZE = 18    # colour-bar title
HEIGHT_FRACTION = 0.3        # x1.0: highest peak = 30 % of smaller span
CONTOUR_STEP_K = 10.0        # contour lines every 10 K (optional)
SMOOTH_FACTOR = 10           # smooth surface: mesh 10x finer than the data
FLOOR_OPACITY = 0.06         # floor: grey at 6 % (nearly black from above)

# Colour bar position (normalized viewport of the 3D view), horizontal
BAR_X, BAR_Y, BAR_W, BAR_H = 0.25, 0.03, 0.50, 0.07

DEFAULT_SIZE = (1600, 900)   # first launch only

# HTML export: Plotly camera vs VTK camera
PLOTLY_FOVY_DEG = 45.0       # Plotly 3D vertical field of view
PLOTLY_SCENE_SCALE = 1.0     # scene units per aspect-ratio unit

PANEL_CSS = """
QScrollArea, QFrame#panel { background: #FAFAF9; }
QScrollArea { border: none; border-left: 1px solid #ddd; }
QLabel.head  { font-weight: bold; font-size: 16px; color: #222; }
QLabel.note  { font-size: 13px; color: #666; }
QLabel.file  { font-size: 14px; color: #085041; }
QLabel.value { font-size: 14px; color: #222; }
QComboBox, QPushButton { font-size: 14px; padding: 6px; }
QCheckBox { font-size: 14px; color: #222; }
"""


# Grid from a HALO .fro: same method B as build_salsa_terrain_grid.py
REST_HZ_DEFAULT = 1420.405752e6
C_KMS = 299792.458
RFI_VMIN_KMS, RFI_VMAX_KMS = 80.0, 95.0  # RFI window (topocentric velocity)
DC_MARGIN_CHANNELS = 2
BASELINE_EDGE_FRACTION = 0.07
GRID_STEP_DEG = 5.0


def peak_method_b(spectrum, velocity):
    """Peak of one spectrum: DC + RFI masked, edge baseline subtracted,
    negative values clipped (as build_salsa_terrain_grid.py, peak_b)."""
    n = len(spectrum)
    spec = spectrum.astype(float)
    mask = np.ones(n, dtype=bool)
    c = n // 2
    mask[max(0, c - DC_MARGIN_CHANNELS):min(n, c + DC_MARGIN_CHANNELS + 1)] = False
    mask[(velocity >= RFI_VMIN_KMS) & (velocity <= RFI_VMAX_KMS)] = False
    edge = int(n * BASELINE_EDGE_FRACTION)
    base = np.zeros(n, dtype=bool)
    base[:edge] = True
    base[-edge:] = True
    base &= mask
    spec_b = np.where(mask, spec - np.mean(spec[base]), np.nan)
    return float(np.nanmax(np.clip(spec_b, 0, None)))


def load_fro_grid(path):
    """GLon x GLat peak grid built on the fly from a HALO .fro file
    (one cell every GRID_STEP_DEG; several spectra in a cell: the highest
    peak, as build_salsa_terrain_grid.py; NaN where nothing observed)."""
    import h5py
    with h5py.File(path, 'r') as f:
        spectra = f['Spectra/raw_spectrum'][:]
        freq = f['Spectra/frequency_axis_hz'][:].astype(float)
        glon = f['Pointing/glon_deg'][:].astype(float)
        glat = f['Pointing/glat_deg'][:].astype(float)
        rest = float(f['Observation'].attrs.get('rest_frequency_hz', 0) or 0)
    if not rest > 0:
        rest = REST_HZ_DEFAULT
    # Topocentric velocity: only used to mask the RFI, which sits at a
    # fixed frequency (same as build_salsa_terrain_grid.py)
    velocity = C_KMS * (rest - freq) / rest
    ok = np.isfinite(glon) & np.isfinite(glat)
    gl = np.round(glon[ok] / GRID_STEP_DEG) * GRID_STEP_DEG
    gb = np.round(glat[ok] / GRID_STEP_DEG) * GRID_STEP_DEG
    peaks = [peak_method_b(s, velocity) for s in spectra[ok]]
    lon_v, lat_v = np.unique(gl), np.unique(gb)
    grid = np.full((len(lat_v), len(lon_v)), np.nan)
    for lo, la, pk in zip(gl, gb, peaks):
        j = int(np.where(lon_v == lo)[0][0])
        i = int(np.where(lat_v == la)[0][0])
        if np.isnan(grid[i, j]) or pk > grid[i, j]:
            grid[i, j] = pk
    print(f'Grid built from {os.path.basename(path)}: {int(ok.sum())} spectra, '
          f'{len(lon_v)} GLon x {len(lat_v)} GLat')
    return lon_v, lat_v, grid


def load_grid(path):
    """Read the grid (.npz from build_salsa_terrain_grid.py, or built from
    a HALO .fro); GLon and GLat sorted in increasing order."""
    if path.lower().endswith('.fro'):
        glon, glat, peak = load_fro_grid(path)
    else:
        d = np.load(path)
        glon = d['glon_values'].astype(float)    # real GLon degrees, 1D
        glat = d['glat_values'].astype(float)    # real GLat degrees, 1D
        peak = d['peak_grid_b'].astype(float)    # (n_glat, n_glon)
    jl, il = np.argsort(glon), np.argsort(glat)
    return glon[jl], glat[il], peak[np.ix_(il, jl)]


def on10(v):
    """True for multiples of 10 deg (graticule)."""
    return abs(v / 10.0 - round(v / 10.0)) < 1e-6


def colorbar_ticks(vmax):
    """Round ticks every 'step' (about four intervals; step 1, 2, 2.5, 3
    or 5 x 10^n, not larger than a quarter of vmax) plus the real maximum;
    a round tick too close to the maximum is dropped (e.g. 0 30 60 90 121)."""
    raw = vmax / 4.0
    mag = 10.0 ** np.floor(np.log10(raw))
    step = max(m * mag for m in (1, 2, 2.5, 3, 5) if m * mag <= raw)
    tk = [float(t) for t in np.arange(0.0, vmax, step)]
    if vmax - tk[-1] < 0.5 * step:
        tk = tk[:-1]
    return tk + [float(round(vmax))]


def smooth_surface(glon, glat, peak, factor, vmax):
    """Cubic spline through the data nodes (exact at the nodes) on a mesh
    'factor' times finer, for display only; clipped to 0..vmax so it
    never shows values beyond the measured range."""
    kx, ky = min(3, len(glat) - 1), min(3, len(glon) - 1)
    spl = RectBivariateSpline(glat, glon, peak, kx=kx, ky=ky, s=0)
    lon_f = np.linspace(glon.min(), glon.max(), (len(glon) - 1) * factor + 1)
    lat_f = np.linspace(glat.min(), glat.max(), (len(glat) - 1) * factor + 1)
    return lon_f, lat_f, np.clip(spl(lat_f, lon_f), 0.0, vmax)


def triangle_mesh(lon, lat, Z, values):
    """Triangle mesh of a (n_lat, n_lon) grid, points row by row. The same
    mesh carries the surface and its contours; graticule and ridge follow
    its edges, so all lines lie exactly on the rendered triangles."""
    nlat, nlon = Z.shape
    X, Y = np.meshgrid(lon, lat)
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    i, j = np.meshgrid(np.arange(nlat - 1), np.arange(nlon - 1),
                       indexing='ij')
    a = (i * nlon + j).ravel()
    b, c, d = a + 1, a + nlon + 1, a + nlon
    three = np.full(a.size, 3)
    faces = np.column_stack([three, a, b, c, three, a, c, d]).ravel()
    mesh = pv.PolyData(pts, faces=faces)
    mesh['Peak T_A'] = values.ravel()
    return mesh


class SalsaTerrainViewer(QMainWindow):

    def __init__(self, path):
        super().__init__()
        self.path = path
        self.setWindowTitle('SALSA survey - peak terrain - HALO')
        glon, glat, peak = load_grid(path)
        print(f'Grid shape: {peak.shape[0]} GLat x {peak.shape[1]} GLon')
        self.n_gaps = int(np.sum(np.isnan(peak)))
        lon_span = float(glon.max() - glon.min())
        lat_span = float(glat.max() - glat.min())
        self.ref_axis, self.ref_span = (('GLat', lat_span)
                                        if lat_span <= lon_span
                                        else ('GLon', lon_span))
        self.peak_max = float(np.nanmax(peak))
        # Display options: all optional
        self.opt = dict(smooth=True, ridge=True, grat=True, cont=False,
                        floor=True, axes=True)

        central = QWidget()
        lay = QHBoxLayout(central)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        self.pl = QtInteractor(central)
        lay.addWidget(self.pl.interactor, 1)
        lay.addWidget(self._build_panel(os.path.basename(path)))
        self.setCentralWidget(central)

        self._build_scene(glon, glat, peak)

        self._settings = QSettings('HALO', 'view_salsa_terrain')
        geom = self._settings.value('geometry')
        if geom is None or not self.restoreGeometry(geom):
            self.resize(*DEFAULT_SIZE)

    # ---------------------------------------------------------------- scene
    def _build_scene(self, glon, glat, peak):
        pl = self.pl
        pl.set_background('black')
        pl.enable_anti_aliasing('ssaa')       # smooth lines, also in the PNG

        # NaN cells are known coverage gaps, not zero emission: shown flat
        # at zero for display only
        self.glon, self.glat = glon, glat
        self.peak = np.nan_to_num(peak, nan=0.0)
        self.vmax = float(self.peak.max())
        lon_span = float(glon.max() - glon.min())
        lat_span = float(glat.max() - glat.min())
        self.base = min(lon_span, lat_span) * HEIGHT_FRACTION
        self.lift = 0.01 * self.base           # lines just above the surface
        self.cbar_ticks = colorbar_ticks(self.vmax)

        # Two versions of surface, ridge and graticule: data nodes ('raw')
        # and cubic spline ('smooth'); the panel shows one of them
        lon_f, lat_f, peak_f = smooth_surface(glon, glat, self.peak,
                                              SMOOTH_FACTOR, self.vmax)
        self.modes = {
            'raw': self._build_mode(glon, glat, self.peak),
            'smooth': self._build_mode(lon_f, lat_f, peak_f),
        }

        self._build_floor(glon, glat)
        self._build_colorbar()
        from vtkmodules.vtkRenderingCore import vtkTextActor
        t = vtkTextActor()
        t.SetInput('SALSA survey - masked + baseline subtracted')
        t.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        t.SetPosition(0.005, 0.99)
        tp = t.GetTextProperty()
        tp.SetFontSize(TITLE_FONT_SIZE)
        tp.SetColor(1.0, 1.0, 1.0)
        tp.SetJustificationToLeft()
        tp.SetVerticalJustificationToTop()
        pl.add_actor(t, reset_camera=False)
        self.title_actor = t
        self._apply_visibility()

        # Initial view: from above, slightly zoomed
        pl.view_xy()
        pl.camera.zoom(1.3)
        self.initial_view = pl.camera_position

    def _build_mode(self, lon, lat, val):
        pl = self.pl
        z0 = val / self.vmax * self.base                  # height x1.0
        mesh = triangle_mesh(lon, lat, z0, val)
        terrain = pl.add_mesh(mesh, scalars='Peak T_A', cmap=BLUE_WHITE,
                              clim=[0.0, self.vmax], smooth_shading=True,
                              show_scalar_bar=False)

        # Galactic plane ridge (b = 0), along a mesh row
        i0 = int(np.argmin(np.abs(lat)))
        ridge_z = z0[i0, :].copy()
        ridge = pv.PolyData(
            np.column_stack([lon, np.zeros(len(lon)), ridge_z + self.lift]),
            lines=np.hstack([[len(lon)], np.arange(len(lon))]))
        ridge_act = pl.add_mesh(ridge, color='yellow', line_width=4,
                                render_lines_as_tubes=True)

        # Graticule on the surface every 10 deg (b = 0 excluded: it is the
        # ridge), along mesh rows and columns
        pts, lines, zb = [], [], []

        def add_line(xs, ys, zs):
            n0 = len(pts)
            for x, y, z in zip(xs, ys, zs):
                pts.append([float(x), float(y), 0.0])
                zb.append(float(z))
            lines.append(np.r_[len(xs), np.arange(n0, n0 + len(xs))])

        for j, lo in enumerate(lon):
            if on10(lo):
                add_line(np.full(len(lat), lo), lat, z0[:, j])
        for i, b in enumerate(lat):
            if on10(b) and abs(b) > 1e-6:
                add_line(lon, np.full(len(lon), b), z0[i, :])
        grat = grat_act = grat_z = None
        if pts:
            grat_z = np.array(zb)
            gp = np.array(pts)
            gp[:, 2] = grat_z + self.lift
            grat = pv.PolyData(gp, lines=np.hstack(lines))
            grat_act = pl.add_mesh(grat, color='white', line_width=1.5,
                                   opacity=0.6)

        return dict(lon=lon, lat=lat, val=val, z0=z0, mesh=mesh,
                    terrain=terrain, ridge=ridge, ridge_z=ridge_z,
                    ridge_act=ridge_act, grat=grat, grat_z=grat_z,
                    grat_act=grat_act)

    def _build_floor(self, glon, glat):
        pl = self.pl
        base = self.base
        m = 0.08 * float(glat.max() - glat.min())    # margin on all sides
        zf = -0.1 * base
        x0, x1 = float(glon.min()) - m, float(glon.max()) + m
        y0, y1 = float(glat.min()) - m, float(glat.max()) + m
        self.floor_geom = (x0, x1, y0, y1, zf, m)

        # Floor plane, slightly larger than the data
        plane = pv.Plane(center=(0.5 * (x0 + x1), 0.5 * (y0 + y1), zf),
                         direction=(0, 0, 1), i_size=x1 - x0,
                         j_size=y1 - y0)
        pl.add_mesh(plane, color='gray', opacity=FLOOR_OPACITY)

        # Floor graticule every 10 deg (white), equator in yellow
        self.lon_ticks = [float(v) for v in
                          np.arange(np.ceil(glon.min() / 10) * 10,
                                    glon.max() + 0.01, 10)]
        self.lat_ticks = [float(v) for v in
                          np.arange(np.ceil(glat.min() / 10) * 10,
                                    glat.max() + 0.01, 10)]
        segs = ([(lo, y0, lo, y1) for lo in self.lon_ticks]
                + [(x0, b, x1, b) for b in self.lat_ticks if abs(b) > 1e-6])
        fpts = []
        for xa, ya, xb, yb in segs:
            fpts += [[xa, ya, zf], [xb, yb, zf]]
        floor_lines = pv.PolyData(
            np.array(fpts, dtype=float),
            lines=np.hstack([[2, 2 * k, 2 * k + 1] for k in range(len(segs))]))
        self.floor_actors = [pl.add_mesh(floor_lines, color='white',
                                         line_width=1, opacity=0.5)]
        if glat.min() <= 0.0 <= glat.max():
            eq = pv.PolyData(np.array([[x0, 0.0, zf], [x1, 0.0, zf]]),
                             lines=np.array([2, 0, 1]))
            self.floor_actors.append(pl.add_mesh(eq, color='yellow',
                                                 line_width=2, opacity=0.8))

        # Axis labels just outside the four floor edges: GLon every 20 deg
        # (near and far), GLat every 10 deg (left and right), plus titles.
        # Billboard text: always faces the camera, drawn from the first
        # frame and in saved PNGs; the terrain can still hide it.
        z_lab = zf + 0.02 * base
        off = 0.6 * m
        lab_pts, lab_txt = [], []
        for lo in np.arange(np.ceil(glon.min() / 20) * 20,
                            glon.max() + 0.01, 20):
            lab_pts += [[lo, y0 - off, z_lab], [lo, y1 + off, z_lab]]
            lab_txt += [f'{lo:.0f}°'] * 2
        for b in self.lat_ticks:
            lab_pts += [[x0 - off, b, z_lab], [x1 + off, b, z_lab]]
            lab_txt += ['0°' if abs(b) < 1e-6 else f'{b:+.0f}°'] * 2
        xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        lab_pts += [[xm, y0 - 2.2 * m, z_lab], [xm, y1 + 2.2 * m, z_lab],
                    [x0 - 3.0 * m, ym, z_lab], [x1 + 3.0 * m, ym, z_lab]]
        lab_txt += (['Galactic longitude l (deg)'] * 2
                    + ['Galactic latitude b (deg)'] * 2)
        from vtkmodules.vtkRenderingCore import vtkBillboardTextActor3D
        self.axes_actors = []
        for (px, py, pz), txt in zip(lab_pts, lab_txt):
            a = vtkBillboardTextActor3D()
            a.SetInput(txt)
            a.SetPosition(px, py, pz)
            tp = a.GetTextProperty()
            tp.SetFontSize(LABEL_FONT_SIZE)
            tp.SetColor(1.0, 1.0, 1.0)
            tp.SetJustificationToCentered()
            tp.SetVerticalJustificationToCentered()
            if txt.startswith('Galactic latitude'):
                tp.SetOrientation(90)            # vertical axis title
            pl.add_actor(a, reset_camera=False)
            self.axes_actors.append(a)

    def _label_font(self, k):
        # 3-D axis labels (pixels): x k while saving a x k PNG, k = 1
        # restores them. The 2-D texts (title, colour-bar title) keep
        # their size: the x k screenshot raises the DPI, which VTK applies
        # only when a text is re-rendered, so they are marked as modified.
        for a in self.axes_actors:
            a.GetTextProperty().SetFontSize(LABEL_FONT_SIZE * k)
        for t, fs in ((self.title_actor, TITLE_FONT_SIZE),
                      (self.cbar_title, CBAR_TITLE_FONT_SIZE)):
            t.GetTextProperty().SetFontSize(fs)
            t.GetTextProperty().Modified()
            t.Modified()

    def _build_colorbar(self):
        # Horizontal colour bar at the bottom, round ticks plus the real
        # maximum. Its title is a separate text centred on the coloured
        # bar (VTK centres its own title on bar + labels).
        pl = self.pl
        sb = pl.add_scalar_bar(title='Peak T_A',
                               mapper=self.modes['raw']['terrain'].mapper,
                               vertical=False, color='white', fmt='%.0f',
                               position_x=BAR_X, position_y=BAR_Y,
                               width=BAR_W, height=BAR_H)
        sb.SetTitle('')
        ticks = vtkDoubleArray()
        for t in self.cbar_ticks:
            ticks.InsertNextValue(t)
        sb.SetUseCustomLabels(True)
        sb.SetCustomLabels(ticks)
        title = vtkTextActor()
        title.SetInput('Peak T_A (K)')
        self.cbar_title = title
        tp = title.GetTextProperty()
        tp.SetColor(1.0, 1.0, 1.0)
        tp.SetFontSize(CBAR_TITLE_FONT_SIZE)
        tp.SetJustificationToCentered()
        tp.SetVerticalJustificationToBottom()
        title.GetPositionCoordinate().SetCoordinateSystemToNormalizedViewport()
        title.SetPosition(BAR_X + BAR_W / 2.0, BAR_Y + BAR_H + 0.01)
        pl.add_actor(title)

    # ------------------------------------------------------------- display
    def _mode(self):
        return self.modes['smooth' if self.opt['smooth'] else 'raw']

    def _scale(self):
        return self.height.value() * HEIGHT_STEP

    def _apply_visibility(self):
        active = 'smooth' if self.opt['smooth'] else 'raw'
        for name, md in self.modes.items():
            on = name == active
            md['terrain'].SetVisibility(on)
            md['ridge_act'].SetVisibility(on and self.opt['ridge'])
            if md['grat_act'] is not None:
                md['grat_act'].SetVisibility(on and self.opt['grat'])
        for act in self.floor_actors:
            act.SetVisibility(self.opt['floor'])
        for a in self.axes_actors:
            a.SetVisibility(self.opt['axes'])
        self._update_contours()
        self.pl.render()

    def _update_contours(self):
        # Contours every CONTOUR_STEP_K kelvin on the shown mesh (height is
        # proportional to Peak T_A), just above the surface
        self.pl.remove_actor('contours')
        if not self.opt['cont']:
            return
        levels = np.arange(CONTOUR_STEP_K, self.vmax, CONTOUR_STEP_K)
        if len(levels) == 0:
            return
        c = self._mode()['mesh'].contour(isosurfaces=levels,
                                         scalars='Peak T_A')
        if c.n_points == 0:
            return
        c.points[:, 2] += self.lift
        self.pl.add_mesh(c, color='#101010', line_width=1, name='contours',
                         reset_camera=False)

    # ---------------------------------------------------------------- panel
    def _build_panel(self, filename):
        panel = QFrame()
        panel.setObjectName('panel')
        box = QVBoxLayout(panel)
        box.setContentsMargins(16, 16, 16, 16)
        box.setSpacing(6)

        def label(text, cls):
            lab = QLabel(text)
            lab.setProperty('class', cls)
            lab.setWordWrap(True)
            return lab

        box.addWidget(label('File', 'head'))
        box.addWidget(label(filename, 'file'))
        box.addSpacing(14)

        # Height scale
        box.addWidget(label('Height scale', 'head'))
        row = QHBoxLayout()
        self.height = QSlider(Qt.Horizontal)
        self.height.setRange(int(round(HEIGHT_MIN / HEIGHT_STEP)),
                             int(round(HEIGHT_MAX / HEIGHT_STEP)))
        self.height.setValue(int(round(1.0 / HEIGHT_STEP)))
        self.height.valueChanged.connect(self._on_height)
        self.height_lab = QLabel('×1.0')
        self.height_lab.setProperty('class', 'value')
        self.height_lab.setFixedWidth(56)
        row.addWidget(self.height)
        row.addWidget(self.height_lab)
        box.addLayout(row)
        box.addWidget(label('Height is relative and arbitrary: the real '
                            'value is the colour.', 'note'))
        box.addWidget(label(f'×1.0: highest peak ({self.peak_max:.0f} K) '
                            f'rises to {HEIGHT_FRACTION * 100:.0f} % of the '
                            f'{self.ref_axis} span ({self.ref_span:.0f}°).',
                            'note'))
        box.addSpacing(14)

        # Display options (all optional)
        box.addWidget(label('Show', 'head'))
        for key, text in (
                ('smooth', 'Smooth surface (spline)'),
                ('ridge', 'Galactic plane (b = 0)'),
                ('grat', 'Graticule on surface'),
                ('cont', f'Contours (every {CONTOUR_STEP_K:.0f} K)'),
                ('floor', 'Graticule on floor'),
                ('axes', 'Axis labels')):
            cb = QCheckBox(text)
            cb.setChecked(self.opt[key])
            cb.toggled.connect(lambda on, k=key: self._on_option(k, on))
            box.addWidget(cb)
        box.addWidget(label('Smooth surface: cubic spline through the data '
                            'nodes, for display only (values at the nodes '
                            'unchanged). Unticked: plain grid of the data '
                            'nodes, linear contours.', 'note'))
        box.addSpacing(14)

        # View
        center = QPushButton('Center view')
        center.clicked.connect(self._on_center)
        box.addWidget(center)
        box.addWidget(label('Frames the terrain, keeping the current '
                            'orientation and height scale.', 'note'))
        box.addSpacing(8)
        reset = QPushButton('Reset view')
        reset.clicked.connect(self._on_reset)
        box.addWidget(reset)
        box.addWidget(label('Restores the initial camera and height '
                            'scale ×1.0.', 'note'))
        box.addSpacing(14)

        # Save PNG
        box.addWidget(label('Save image (PNG)', 'head'))
        row = QHBoxLayout()
        self.res_combo = QComboBox()
        for k in (1, 2, 3):
            self.res_combo.addItem(f'×{k} resolution', k)
        self.res_combo.setCurrentIndex(0)                 # x1 by default
        save = QPushButton('Save PNG')
        save.clicked.connect(self._on_save)
        row.addWidget(self.res_combo)
        row.addWidget(save)
        box.addLayout(row)
        self.status = label('Opens a "Save as" window in the dataset '
                            'folder; nothing is written if you cancel. ×2 '
                            'and ×3 multiply the window size, for sharper '
                            'zooming.', 'note')
        box.addWidget(self.status)
        box.addSpacing(14)

        # Interactive HTML page
        box.addWidget(label('Interactive page (HTML)', 'head'))
        row = QHBoxLayout()
        prev = QPushButton('Preview in browser')
        prev.clicked.connect(lambda: self._on_html(preview=True))
        html = QPushButton('Export HTML…')
        html.clicked.connect(lambda: self._on_html(preview=False))
        row.addWidget(prev)
        row.addWidget(html)
        box.addLayout(row)
        self.html_status = label('Self-contained page (Plotly) to rotate '
                                 'in any browser, also offline. Opens with '
                                 'the current view, height scale and '
                                 'options. Preview: not saved.', 'note')
        box.addWidget(self.html_status)

        box.addStretch(1)
        box.addSpacing(14)
        box.addWidget(label('Peak T_A: antenna temperature, nominal '
                            '(not calibrated).', 'note'))
        box.addWidget(label('DC spike and RFI masked, per-spectrum baseline '
                            'subtracted (method B).', 'note'))
        if self.n_gaps:
            box.addWidget(label(f'{self.n_gaps} grid cells without '
                                'observations are shown flat at zero '
                                '(coverage gaps, not zero emission).', 'note'))

        # Scroll bar when the window is not tall enough
        scroll = QScrollArea()
        scroll.setWidget(panel)
        scroll.setWidgetResizable(True)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        scroll.setFixedWidth(350)
        scroll.setStyleSheet(PANEL_CSS)
        return scroll

    # ------------------------------------------------------------ callbacks
    def _on_option(self, key, on):
        self.opt[key] = bool(on)
        self._apply_visibility()

    def _on_height(self, value):
        s = value * HEIGHT_STEP
        self.height_lab.setText(f'×{s:.1f}')
        for md in self.modes.values():
            md['mesh'].points[:, 2] = (md['z0'] * s).ravel()
            md['ridge'].points[:, 2] = md['ridge_z'] * s + self.lift
            if md['grat'] is not None:
                md['grat'].points[:, 2] = md['grat_z'] * s + self.lift
        self._update_contours()
        self.pl.render()

    def _on_center(self):
        # ResetCamera frames the visible actors, keeping orientation
        self.pl.reset_camera()
        self.pl.render()

    def _on_reset(self):
        self.height.setValue(int(round(1.0 / HEIGHT_STEP)))   # also redraws
        self.pl.camera_position = self.initial_view
        self.pl.render()

    def _ask_save_path(self, ext, file_filter):
        # "Save as" window in the dataset folder with a proposed name
        # (date and time); returns '' if cancelled. Qt asks before
        # overwriting an existing file.
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        start = os.path.join(os.path.dirname(os.path.abspath(self.path)),
                             f'salsa_terrain_{stamp}.{ext}')
        out, _ = QFileDialog.getSaveFileName(
            self, f'Save {ext.upper()}', start, file_filter,
            options=QFileDialog.DontUseNativeDialog)
        if out and not out.lower().endswith('.' + ext):
            out += '.' + ext
        return out

    def _on_save(self):
        k = int(self.res_combo.currentData())
        out = self._ask_save_path('png', 'PNG image (*.png)')
        if not out:
            return
        try:
            self._label_font(k)
            try:
                img = self.pl.screenshot(out, scale=k, return_img=True)
            finally:
                self._label_font(1)
        except Exception as e:
            QMessageBox.warning(self, 'Save failed',
                                f'Could not save the PNG image:\n{e}')
            return
        h, w = img.shape[:2]
        self.status.setText(f'Saved: {os.path.basename(out)} '
                            f'({w} × {h} px, ×{k})')
        print(f'Saved: {out} ({w} x {h} px)')

    def _on_html(self, preview=False):
        # Same surface (smooth or raw), colours, options and proportions as
        # the viewer, at the current height scale and from the current
        # view. Plotly's own 3D axes are not used (Plotly moves them
        # between box edges): floor, graticule and labels are drawn as in
        # the viewer. plotly.js is embedded (offline, shareable).
        import plotly.graph_objects as go
        s = self._scale()
        md = self._mode()
        lon, lat = md['lon'], md['lat']
        x0, x1, y0, y1, zf, m = self.floor_geom
        flat = FLAT_LIGHT

        # Surface; fixed contours only if ticked; no lines following the
        # mouse; flat lighting (the colour is the data)
        k2z = self.base * s / self.vmax
        no_hl = dict(highlight=False)
        contours = dict(x=no_hl, y=no_hl, z=no_hl)
        if self.opt['cont']:
            contours = dict(x=no_hl, y=no_hl,
                            z=dict(show=True, highlight=False,
                                   start=CONTOUR_STEP_K * k2z,
                                   end=self.vmax * k2z,
                                   size=CONTOUR_STEP_K * k2z,
                                   color='#101010', width=1,
                                   usecolormap=False))
        X, Y = np.meshgrid(lon, lat)
        traces = [go.Surface(
            x=X, y=Y, z=md['z0'] * s, surfacecolor=md['val'],
            cmin=0.0, cmax=self.vmax, colorscale=PLOTLY_COLORSCALE,
            lighting=flat, contours=contours, hoverinfo='skip',
            colorbar=dict(title=dict(text='Peak T_A (K)', side='top'),
                          orientation='h', x=0.5, xanchor='center',
                          y=0.02, yanchor='bottom', len=0.5, thickness=18,
                          tickformat='.0f', tickvals=self.cbar_ticks))]

        # Hover values: invisible markers on the data nodes only (measured
        # values; the surface does not pass its values to the hover box)
        Xn, Yn = np.meshgrid(self.glon, self.glat)
        traces.append(go.Scatter3d(
            x=Xn.ravel(), y=Yn.ravel(),
            z=(self.peak / self.vmax * self.base * s).ravel(),
            mode='markers', marker=dict(size=6, opacity=0.01, color='white'),
            customdata=self.peak.ravel(),
            hovertemplate=('l = %{x:.1f}°<br>b = %{y:.1f}°<br>'
                           'Peak T_A = %{customdata:.1f} K<extra></extra>')))

        # Floor (translucent, as in the viewer)
        traces.append(go.Surface(
            x=[[x0, x1], [x0, x1]], y=[[y0, y0], [y1, y1]],
            z=[[zf, zf], [zf, zf]], surfacecolor=[[0, 0], [0, 0]],
            colorscale=[[0, '#808080'], [1, '#808080']], showscale=False,
            opacity=FLOOR_OPACITY, lighting=flat, hoverinfo='skip'))

        # Floor graticule: opaque grey (translucent lines ignore depth in
        # Plotly and would show through the terrain), equator in yellow
        if self.opt['floor']:
            fx, fy, fz = [], [], []
            for lo in self.lon_ticks:
                fx += [lo, lo, None]; fy += [y0, y1, None]
                fz += [zf, zf, None]
            for b in self.lat_ticks:
                if abs(b) > 1e-6:
                    fx += [x0, x1, None]; fy += [b, b, None]
                    fz += [zf, zf, None]
            traces.append(go.Scatter3d(x=fx, y=fy, z=fz, mode='lines',
                                       hoverinfo='skip',
                                       line=dict(color='#808080', width=1)))
            if self.glat.min() <= 0.0 <= self.glat.max():
                traces.append(go.Scatter3d(
                    x=[x0, x1], y=[0.0, 0.0], z=[zf, zf], mode='lines',
                    hoverinfo='skip', line=dict(color='yellow', width=3)))

        # Axis labels on the four floor edges
        if self.opt['axes']:
            lab_lon = [lo for lo in self.lon_ticks
                       if abs(lo / 20 - round(lo / 20)) < 1e-6]
            nl, nb = len(lab_lon), len(self.lat_ticks)
            tx = lab_lon * 2 + [x0 - 0.6 * m] * nb + [x1 + 0.6 * m] * nb
            ty = ([y0 - 0.6 * m] * nl + [y1 + 0.6 * m] * nl
                  + self.lat_ticks * 2)
            tt = ([f'{lo:.0f}°' for lo in lab_lon] * 2
                  + ['0°' if abs(b) < 1e-6 else f'{b:+.0f}°'
                     for b in self.lat_ticks] * 2)
            xm, ym = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
            tx += [xm, xm, x0 - 3.0 * m, x1 + 3.0 * m]
            ty += [y0 - 2.2 * m, y1 + 2.2 * m, ym, ym]
            tt += (['Galactic longitude l (deg)'] * 2
                   + ['Galactic latitude b (deg)'] * 2)
            traces.append(go.Scatter3d(x=tx, y=ty, z=[zf] * len(tx),
                                       mode='text', text=tt,
                                       hoverinfo='skip',
                                       textfont=dict(color='white', size=12)))

        # Galactic plane ridge and surface graticule, as shown
        if self.opt['ridge']:
            traces.append(go.Scatter3d(
                x=lon, y=np.zeros_like(lon), z=md['ridge_z'] * s + self.lift,
                mode='lines', line=dict(color='yellow', width=6),
                hoverinfo='skip'))
        if self.opt['grat'] and md['grat'] is not None:
            gx, gy, gz = [], [], []
            for j, lo in enumerate(lon):
                if on10(lo):
                    gx += [float(lo)] * len(lat) + [None]
                    gy += [float(b) for b in lat] + [None]
                    gz += [float(z) for z in
                           md['z0'][:, j] * s + self.lift] + [None]
            for i, b in enumerate(lat):
                if on10(b) and abs(b) > 1e-6:
                    gx += [float(v) for v in lon] + [None]
                    gy += [float(b)] * len(lon) + [None]
                    gz += [float(z) for z in
                           md['z0'][i, :] * s + self.lift] + [None]
            traces.append(go.Scatter3d(
                x=gx, y=gy, z=gz, mode='lines', hoverinfo='skip',
                line=dict(color='#d8d8d8', width=2)))     # opaque (depth)

        # Camera: the page opens with the viewer's current view (what you
        # see is what you get). Explicit axis ranges and a data-proportional
        # aspect ratio make Plotly's scene coordinates a uniform scaling of
        # the data coordinates, so the VTK camera converts exactly; the
        # distance is corrected for the different field of view.
        ztop = float((md['z0'] * s).max()) + 0.05 * self.base
        rng = [(x0 - 3.5 * m, x1 + 3.5 * m), (y0 - 2.5 * m, y1 + 2.5 * m),
               (zf, max(ztop, zf + 1e-6))]
        L = np.array([r[1] - r[0] for r in rng])
        mid = np.array([0.5 * (r[0] + r[1]) for r in rng])
        ar = L / L[0]
        k = PLOTLY_SCENE_SCALE / L[0]          # data units -> scene units
        cam = self.pl.camera
        pos = np.array(cam.position, dtype=float)
        foc = np.array(cam.focal_point, dtype=float)
        upv = np.array(cam.up, dtype=float)
        fov = (np.tan(np.radians(cam.view_angle) / 2.0)
               / np.tan(np.radians(PLOTLY_FOVY_DEG) / 2.0))
        c_n = (foc - mid) * k
        e_n = c_n + (pos - foc) * k * fov

        def xyz(v):
            return dict(x=float(v[0]), y=float(v[1]), z=float(v[2]))

        fig = go.Figure(traces)
        fig.update_layout(
            title=dict(text=('SALSA survey - peak terrain '
                             '(masked + baseline subtracted)<br><sup>'
                             'Peak T_A: antenna temperature, nominal (not '
                             'calibrated). Height is relative and arbitrary: '
                             'the real value is the colour.</sup>'),
                       font=dict(color='white')),
            paper_bgcolor='black', font=dict(color='white'),
            showlegend=False, margin=dict(l=0, r=0, t=80, b=10),
            scene=dict(xaxis=dict(visible=False, range=list(rng[0])),
                       yaxis=dict(visible=False, range=list(rng[1])),
                       zaxis=dict(visible=False, range=list(rng[2])),
                       aspectmode='manual',
                       aspectratio=dict(x=float(ar[0]), y=float(ar[1]),
                                        z=float(ar[2])),
                       camera=dict(eye=xyz(e_n), center=xyz(c_n),
                                   up=xyz(upv)),
                       bgcolor='black',
                       domain=dict(x=[0, 1], y=[0.14, 1])))

        if preview:
            # Not saved in the project: one temporary file in the system
            # temp folder, overwritten at every preview
            out = os.path.join(tempfile.gettempdir(),
                               'halo_salsa_terrain_preview.html')
        else:
            out = self._ask_save_path('html', 'HTML page (*.html)')
            if not out:
                return
        try:
            fig.write_html(out, include_plotlyjs=True, full_html=True,
                           config={'displaylogo': False,
                                   'toImageButtonOptions': {'format': 'png',
                                                            'scale': 3}})
        except Exception as e:
            QMessageBox.warning(self, 'Export failed',
                                f'Could not write the HTML page:\n{e}')
            return
        what = ('Preview (not saved)' if preview
                else f'Saved: {os.path.basename(out)}')
        self.html_status.setText(f'{what} (height ×{s:.1f}), '
                                 f'opened in the browser')
        print(f'HTML: {out}')
        QDesktopServices.openUrl(QUrl.fromLocalFile(out))

    def closeEvent(self, event):
        self._settings.setValue('geometry', self.saveGeometry())
        self.pl.close()
        super().closeEvent(event)


def main():
    if len(sys.argv) < 2:
        sys.exit('Usage: view_salsa_terrain.py <file.fro | salsa_terrain_grid.npz>')
    path = sys.argv[1]
    if not os.path.isfile(path):
        sys.exit(f'File not found: {path}')
    app = QApplication(sys.argv)
    win = SalsaTerrainViewer(path)
    win.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
