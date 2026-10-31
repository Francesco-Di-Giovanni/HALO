#!/usr/bin/env python3
"""
fro_onoff_stack.py - FRO ON/OFF Stack (OOS), HALO module.

Position-switched spectral line from one or more .fro files, one file per
observing session (FRO v1.8 or later: the files must carry the per-spectrum
LSR correction Doppler/vlsr_correction_mps and the rest frequency of the
observed line Observation/rest_frequency_hz).

Method, for every ON/OFF pair:  T_A = Tsys * (ON - OFF) / OFF, on an LSR
velocity axis for the observed line, using the LSR correction of that
pair's ON spectrum only. The pairs are aligned on a common LSR axis by
whole-channel shifts (no interpolation, so the noise statistics are
untouched; residual misalignment <= half a channel) and summed with
radiometer weights. A polynomial baseline is fitted outside the expected
line window and the measured noise is compared with the radiometer
equation.

Usage:
  fro_onoff_stack.py [options] <file.fro> [<file.fro> ...]
  fro_onoff_stack.py [options] --search-dir DIR [--prefer FILE]
  Without files and without --search-dir a dialog opens (multiple
  selection inside one folder).

Options:
  --search-dir DIR  list every .fro under DIR and choose the sessions from a
                    numbered list (this is how the HALO Launcher starts it)
  --prefer FILE     with --search-dir: the proposed selection is every file
                    with the same target and rest frequency as FILE
  --channel-axis    channel axis on top of both panels and highlighted
                    watched spur channels (instrument table)
  --inject MK       TEST ONLY: add synthetic peaks of MK mK at the velocities
                    listed for the target in known_targets.json
  --yes             accept every proposed answer, no questions
  --no-open         do not open the SVG at the end (it opens by default)

Questions (type < and Enter to go back one step): sessions to stack (with
--search-dir), which pointing is the OFF position, the expected line
window, the nominal Tsys (only for facilities not in the instrument table),
the output folder.

Known targets (line window, peak velocities, reference) are read from
known_targets.json next to this file; instrument parameters (Tsys, masked
edge and band-centre channels, watched spur channels) from the INSTRUMENTS
table below.

Outputs: PNG, SVG (with every channel number along the bottom of each
panel, readable when zoomed in) and NPZ with the summed spectrum.

Channels: 0-based. The "reference grid" is the channel grid of the sum; a
channel c of a pair lands on reference channel c - shift.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import textwrap
import warnings

import h5py
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

warnings.filterwarnings('ignore')

C = 299792.458              # km/s
HERE = os.path.dirname(os.path.abspath(__file__))
TARGETS_FILE = os.path.join(HERE, 'known_targets.json')
BACK = '<'
AUTO = False                # set by --yes

# Instrument parameters, by Source/facility (upper case).
INSTRUMENTS = {
    'SALSA': dict(
        tsys=300.0,         # K, nominal (approximate, see SALSA docs)
        edge=40,            # channels masked at each band edge
        dc_half=1,          # channels masked each side of the band centre
        mask_dc=True,
        spur_mask=(),       # channels always masked (none confirmed yet)
        spur_watch=(320, 322),  # highlighted only, with --channel-axis
        note='SALSA nominal Tsys; band-centre LO leakage masked'),
}
NO_INSTRUMENT = dict(edge=0, dc_half=0, mask_dc=False, spur_mask=(),
                     spur_watch=(), note='facility not in the instrument '
                     'table: no masks')

OFF_TOL = 1.0               # deg, |GLAT - OFF latitude| for OFF spectra
INJ_FWHM = 1.5              # km/s, width of the synthetic test peaks
BASE_DEG = 3                # baseline polynomial degree
RULER_FONT = 2.5            # pt, channel numbers in the SVG ruler
FREQ_AXIS_POS = -0.22       # frequency axis position, axes fraction
FREQ_COLOR = 'tab:blue'     # colour of the frequency axis
NOTE_FONT = 13              # pt, notes below the figure
NOTE_CHAR_EM = 0.56         # average character width (em), to wrap the notes
FIG_SIZE = (21, 12)         # figure size, inches (width, height)
SVG_VIEWERS = ('loupe', 'eog')  # tried in order (fit to window, vector zoom); else xdg-open

plt.rcParams['font.size'] = 12


def dec(x):
    return x.decode() if isinstance(x, bytes) else str(x)


def ask(prompt, default):
    """Ask for a value showing the default; returns BACK for '<'."""
    if AUTO:
        if default == '':
            sys.exit(f'ERROR: --yes given but no proposed answer for: {prompt}')
        print(f'{prompt} [{default}]: (accepted)')
        return default
    val = input(f'{prompt} [{default}]: ').strip()
    if val == BACK:
        return BACK
    return val if val else default


def parse_selection(text, max_n):
    """Printer-style selection ('1-3,5') -> sorted unique indices (1-based).
    Raises ValueError on malformed input."""
    idx = set()
    for part in text.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            a, b = part.split('-', 1)
            idx.update(range(int(a), int(b) + 1))
        else:
            idx.add(int(part))
    return sorted(i for i in idx if 1 <= i <= max_n)


def compress(idx):
    """[1,2,3,5] -> '1-3,5'."""
    out, i = [], 0
    while i < len(idx):
        j = i
        while j + 1 < len(idx) and idx[j + 1] == idx[j] + 1:
            j += 1
        out.append(str(idx[i]) if i == j else f'{idx[i]}-{idx[j]}')
        i = j + 1
    return ','.join(out)


def pick_files():
    """Multiple-selection file dialog (Qt, zenity fallback)."""
    start = os.path.expanduser('~/FRO/Parsers')
    try:
        from PyQt5.QtWidgets import QApplication, QFileDialog
        app = QApplication.instance() or QApplication(sys.argv[:1])
        files, _ = QFileDialog.getOpenFileNames(
            None, 'Select one or more .fro files (one per session)', start,
            'FRO files (*.fro)')
        return list(files)
    except ImportError:
        r = subprocess.run(['zenity', '--file-selection', '--multiple',
                            '--separator=\n', '--file-filter=*.fro',
                            f'--filename={start}/'],
                           capture_output=True, text=True)
        return [x for x in r.stdout.splitlines() if x]


def scan_candidates(root):
    """Light summary of every .fro under root (attributes only)."""
    out = []
    for dp, dn, fn in os.walk(root):
        dn.sort()
        for f in sorted(fn):
            if not f.endswith('.fro'):
                continue
            p = os.path.join(dp, f)
            e = dict(path=p, target='?', rest=np.nan, n=0, start='',
                     version='?', lsr=False)
            try:
                with h5py.File(p, 'r') as h:
                    obs = h['Observation'].attrs
                    e.update(
                        target=dec(obs.get('target', '')),
                        rest=float(obs.get('rest_frequency_hz', np.nan)),
                        n=h['Spectra/raw_spectrum'].shape[0],
                        start=dec(h['Time'].attrs.get('start_utc', ''))[:10],
                        version=dec(h.attrs.get('format_version', '?')),
                        lsr='Doppler/vlsr_correction_mps' in h)
            except Exception:
                pass
            out.append(e)
    out.sort(key=lambda e: (e['start'], e['path']))
    return out


def load_fro(path):
    with h5py.File(path, 'r') as f:
        n = f['Spectra/raw_spectrum'].shape[0]

        def rows(p):
            return (f[p][...].astype(float) if p in f
                    else np.full(n, np.nan))

        obs = f['Observation'].attrs
        s = dict(
            path=path,
            label=os.path.splitext(os.path.basename(path))[0],
            version=dec(f.attrs.get('format_version', '?')),
            data=f['Spectra/raw_spectrum'][...].astype(float),
            freq=f['Spectra/frequency_axis_hz'][...].astype(float),
            glat=rows('Pointing/glat_deg'),
            el=rows('Pointing/elevation_deg'),
            t=rows('Time/integration_time_s'),
            unix=rows('Time/unix_time'),
            utc=[dec(x) for x in f['Time/utc'][...]],
            vsys=(f['Doppler/vlsr_correction_mps'][...].astype(float) / 1000.0
                  if 'Doppler/vlsr_correction_mps' in f else None),
            names=([dec(x) for x in f['Source/file_name'][...]]
                   if 'Source/file_name' in f else [''] * n),
            rest=float(obs.get('rest_frequency_hz', np.nan)),
            target=dec(obs.get('target', '')),
            facility=dec(f['Source'].attrs.get('facility', '')),
            instrument=dec(f['Source'].attrs.get('instrument', '')))
    s['names'] = [nm if nm else f'row {i}' for i, nm in enumerate(s['names'])]
    return s


def check_sessions(sess):
    """Return an error message, or None if the sessions can be stacked."""
    for s in sess:
        tag = os.path.basename(s['path'])
        if s['vsys'] is None or np.all(np.isnan(s['vsys'])):
            return (f'ERROR: {tag} (FRO v{s["version"]}) has no per-spectrum '
                    'LSR correction (Doppler/vlsr_correction_mps). Convert the '
                    'session again with a parser that writes FRO v1.8.')
        if np.any(np.isnan(s['vsys'])):
            return f'ERROR: {tag}: LSR correction missing for some spectra.'
        if not np.isfinite(s['rest']):
            return (f'ERROR: {tag}: no rest frequency of the observed line '
                    '(Observation/rest_frequency_hz).')
        if np.any(~np.isfinite(s['t'])):
            return f'ERROR: {tag}: integration time missing for some spectra.'
    ref = sess[0]
    for s in sess[1:]:
        if (s['freq'].shape != ref['freq'].shape
                or not np.array_equal(s['freq'], ref['freq'])):
            return ('ERROR: frequency axes differ between files (all sessions '
                    'must use the same settings).')
        if s['rest'] != ref['rest']:
            return 'ERROR: the files refer to different rest frequencies.'
    targets = sorted({s['target'] for s in sess})
    if len(targets) > 1:
        print(f'WARNING: different targets in the selection: {targets}')
    return None


def prepare(st):
    """Derive target, instrument and proposed answers from the selected
    sessions. Return an error message, or None."""
    sess = st['sess']
    key = tuple(s['path'] for s in sess)
    if key == st['key']:
        return None
    ref = sess[0]
    allb = np.round(np.concatenate([s['glat'] for s in sess]), 1)
    bvals, bcount = np.unique(allb[np.isfinite(allb)], return_counts=True)
    if bvals.size < 2:
        return ('ERROR: only one pointing found in the selection: ON and OFF '
                'positions are needed.')
    st['key'] = key
    st['ref'] = ref
    st['facility'], st['instrument'] = ref['facility'], ref['instrument']
    st['target'], st['rest'] = ref['target'], ref['rest']
    try:
        known = json.load(open(TARGETS_FILE))
    except (OSError, ValueError):
        known = {}
    tinfo = known.get(st['target'], {})
    st['peaks'] = [tuple(p) for p in tinfo.get('peaks_kms', [])]
    st['peak_ref'] = tinfo.get('reference', '')
    st['inst'] = INSTRUMENTS.get(st['facility'].upper())
    st['bvals'], st['bcount'] = bvals, bcount
    st['off_choice'] = str(int(np.argmax(np.abs(bvals))) + 1)
    win0 = tinfo.get('line_window_kms')
    if not st['win_custom']:
        st['win_str'] = f'{win0[0]:g},{win0[1]:g}' if win0 else ''
    if st['inst']:
        st['tsys_str'] = f'{st["inst"]["tsys"]:g}'
    st['multi'] = len(sess) > 1
    st['tag'] = (st['target'] or 'target').replace(' ', '_')
    if not st['out_custom']:
        if st['multi']:
            common = os.path.commonpath([os.path.dirname(s['path'])
                                         for s in sess])
            st['out'] = os.path.join(common,
                                     f'{st["facility"]}_{st["tag"]}_stack')
        else:
            st['out'] = os.path.dirname(ref['path'])
    print(f'\n{len(sess)} file(s): {st["facility"]} {st["instrument"]}, '
          f'target "{st["target"] or "-"}", rest frequency '
          f'{st["rest"] / 1e6:.6f} MHz')
    print(f'Known target: '
          f'{"yes (" + st["peak_ref"] + ")" if tinfo else "no"}; '
          f'instrument table: {"yes" if st["inst"] else "no (no masks)"}')
    return None


def stamp(name):
    return os.path.splitext(name)[0].split('-')[-1]


def make_pairs(s, b_off):
    """ON/OFF pairs of one session (file), in time order."""
    order = [int(i) for i in np.argsort(s['unix'], kind='stable')]
    is_off = np.abs(s['glat'] - b_off) <= OFF_TOL
    print(f'\nSession {s["label"]}')
    print(f'  {"spectrum":36s} {"GLAT":>6s} {"El":>6s} {"t (s)":>6s} '
          f'{"VLSRcor":>8s}  role')
    for i in order:
        print(f'  {s["names"][i]:36s} {s["glat"][i]:6.1f} {s["el"][i]:6.1f} '
              f'{s["t"][i]:6.0f} {s["vsys"][i]:+8.2f}  '
              f'{"OFF" if is_off[i] else "ON"}')
    on = [i for i in order if not is_off[i]]
    off = [i for i in order if is_off[i]]
    if not on or not off:
        sys.exit(f'ERROR: {s["label"]}: need both ON and OFF '
                 f'({len(on)} ON, {len(off)} OFF)')
    pos = {i: k for k, i in enumerate(order)}
    if len(on) == len(off):
        match = list(zip(on, off))
    else:
        print(f'  WARNING: {len(on)} ON and {len(off)} OFF: each ON paired '
              'with the nearest OFF in time (radiometer estimate optimistic)')
        match = [(a, min(off, key=lambda b: abs(pos[b] - pos[a]))) for a in on]
    return [dict(session=s['label'], on_data=s['data'][a], off_data=s['data'][b],
                 t_on=s['t'][a], t_off=s['t'][b], vsys=s['vsys'][a],
                 on_name=s['names'][a], off_name=s['names'][b])
            for a, b in match]


def channel_mask(n, inst):
    good = np.ones(n, bool)
    if inst['edge']:
        good[:inst['edge']] = good[-inst['edge']:] = False
    if inst['mask_dc']:
        dc = n // 2
        good[dc - inst['dc_half']:dc + inst['dc_half'] + 1] = False
    for c in inst['spur_mask']:
        good[c] = False
    return good


def stack(pairs, n):
    """Weighted mean on the reference channel grid (NaN where no pair)."""
    num = np.zeros(n)
    wsum = np.zeros(n)
    j = np.arange(n)
    for p in pairs:
        i = j + p['shift']
        ok = (i >= 0) & (i < n)
        ta = np.full(n, np.nan)
        ta[ok] = p['ta'][i[ok]]
        use = np.isfinite(ta)
        num[use] += p['w'] * ta[use]
        wsum[use] += p['w']
    with np.errstate(invalid='ignore', divide='ignore'):
        return num / wsum, wsum


def baseline(v, ta, win):
    """Polynomial baseline outside the line window; robust noise.
    Narrow interference is searched for only outside the line window."""
    n = ta.size
    valid = np.isfinite(ta)
    line = (v > win[0]) & (v < win[1])
    fit = valid & ~line
    x = (np.arange(n) - n / 2) / (n / 2)
    base = np.polyval(np.polyfit(x[fit], ta[fit], BASE_DEG), x)
    res = ta - base
    mad = 1.4826 * np.median(np.abs(res[fit] - np.median(res[fit])))
    spikes = fit & (np.abs(np.nan_to_num(res)) > 5.0 * mad)
    rms = np.std(res[fit & ~spikes])
    return res, base, line, spikes, rms


# ---------------------------------------------------------------------------
# Interactive steps: each returns BACK to go to the previous one, else None.
# ---------------------------------------------------------------------------

def step_files(st):
    c = st['cands']
    print(f'\n.fro files found under {st["root"]}:')
    for k, e in enumerate(c, 1):
        rel = os.path.relpath(e['path'], st['root'])
        rest = (f'{e["rest"] / 1e6:.3f} MHz' if np.isfinite(e['rest'])
                else 'rest ?')
        flag = '' if e['lsr'] else f'   (no LSR correction, FRO v{e["version"]})'
        print(f'  {k:2d}. {rel}   [{e["target"] or "-"}, {e["start"]}, '
              f'{e["n"]} spectra, {rest}]{flag}')
    while True:
        v = ask("Select the sessions to stack (e.g. '1-5' or '1,3')",
                st['sel_str'])
        if v == BACK:
            return BACK
        try:
            idx = parse_selection(v, len(c))
        except ValueError:
            idx = []
        err = None
        if not idx:
            err = 'Not a valid selection.'
        else:
            sess = [load_fro(c[i - 1]['path']) for i in idx]
            sess.sort(key=lambda s: np.nanmin(s['unix']))
            err = check_sessions(sess)
            if not err:
                st['sess'] = sess
                err = prepare(st)
        if err:
            if AUTO:
                sys.exit(err)
            print(f'  {err} Try again.')
            continue
        st['sel_str'] = v
        return None


def step_off(st):
    print('\nPointings found (galactic latitude):')
    for k, (b, c) in enumerate(zip(st['bvals'], st['bcount']), 1):
        print(f'  {k}. b = {b:+.1f} deg  ({c} spectra)')
    nb = st['bvals'].size
    while True:
        v = ask('Which one is the OFF position? Choose a number',
                st['off_choice'])
        if v == BACK:
            return BACK
        if v.isdigit() and 1 <= int(v) <= nb:
            st['off_choice'] = v
            st['b_off'] = float(st['bvals'][int(v) - 1])
            return None
        if AUTO:
            sys.exit('ERROR: invalid proposed answer.')
        print(f'  "{v}" is not valid: type a number from 1 to {nb}.')


def step_window(st):
    while True:
        v = ask('Expected line window, LSR velocity min,max (km/s)',
                st['win_str'])
        if v == BACK:
            return BACK
        try:
            vals = [float(x) for x in v.replace(';', ',').split(',')]
        except ValueError:
            vals = []
        if len(vals) != 2 or vals[0] == vals[1]:
            if AUTO:
                sys.exit('ERROR: invalid line window.')
            print('  Type two different numbers separated by a comma, '
                  'e.g. -75,-35.')
            continue
        lo, hi = sorted(vals)
        new = f'{lo:g},{hi:g}'
        st['win_custom'] = st['win_custom'] or (new != st['win_str'])
        st['win_str'] = new
        st['win'] = (lo, hi)
        return None


def step_tsys(st):
    if st['inst']:
        st['tsys'] = st['inst']['tsys']
        return None
    while True:
        v = ask('Nominal Tsys (K)', st['tsys_str'])
        if v == BACK:
            return BACK
        try:
            t = float(v)
        except ValueError:
            t = -1
        if t <= 0:
            if AUTO:
                sys.exit('ERROR: invalid Tsys.')
            print('  Type a positive number.')
            continue
        st['tsys_str'] = f'{t:g}'
        st['tsys'] = t
        return None


def step_out(st):
    v = ask('Output folder', st['out'])
    if v == BACK:
        return BACK
    new = os.path.abspath(os.path.expanduser(v))
    st['out_custom'] = st['out_custom'] or (new != st['out'])
    st['out'] = new
    return None


def run_steps(steps, st):
    i = 0
    while i < len(steps):
        if steps[i](st) == BACK:
            if i == 0:
                print('  Already at the first question.')
            else:
                i -= 1
            continue
        i += 1


def main():
    global AUTO
    ap = argparse.ArgumentParser(
        description='FRO ON/OFF Stack: position-switched line from one or '
                    'more .fro sessions.')
    ap.add_argument('files', nargs='*', help='.fro file(s), one per session')
    ap.add_argument('--search-dir', help='choose the sessions among the .fro '
                    'files under this folder')
    ap.add_argument('--prefer', help='with --search-dir: propose the files '
                    'with the same target as this one')
    ap.add_argument('--channel-axis', action='store_true',
                    help='channel axis on top, highlighted watched spur '
                         'channels')
    ap.add_argument('--inject', type=float, default=0.0,
                    help='TEST ONLY: synthetic peak amplitude in mK')
    ap.add_argument('--yes', action='store_true',
                    help='accept every proposed answer, no questions')
    ap.add_argument('--no-open', action='store_true',
                    help='do not open the SVG at the end')
    a = ap.parse_args()
    AUTO = a.yes

    st = dict(sess=None, key=None, cands=[], root='', sel_str='',
              off_choice='', b_off=None, win_str='', win=None,
              win_custom=False, tsys_str='', tsys=None, out='',
              out_custom=False)
    steps = [step_off, step_window, step_tsys, step_out]

    if not a.files and a.search_dir:
        root = os.path.abspath(os.path.expanduser(a.search_dir))
        cands = scan_candidates(root)
        if not cands:
            sys.exit(f'No .fro files under {root}.')
        st['cands'], st['root'] = cands, root
        if a.prefer:
            pref = next((e for e in cands
                         if os.path.abspath(e['path']) == os.path.abspath(a.prefer)),
                        None)
            if pref:
                idx = [k for k, e in enumerate(cands, 1)
                       if e['lsr'] and e['target'] == pref['target']
                       and e['rest'] == pref['rest']]
                st['sel_str'] = compress(idx) if idx else ''
        steps = [step_files] + steps
    else:
        files = a.files or pick_files()
        if not files:
            sys.exit('No .fro file selected.')
        files = [os.path.abspath(f) for f in files]
        for f in files:
            if not os.path.isfile(f):
                sys.exit(f'ERROR: not a file: {f}')
        sess = [load_fro(f) for f in files]
        sess.sort(key=lambda s: np.nanmin(s['unix']))
        err = check_sessions(sess)
        if not err:
            st['sess'] = sess
            err = prepare(st)
        if err:
            sys.exit(err)

    if not AUTO:
        print('\nType < and Enter at any question to go back one step.')
    run_steps(steps, st)

    sess, ref = st['sess'], st['ref']
    facility, instrument = st['facility'], st['instrument']
    target, rest, peaks, peak_ref = (st['target'], st['rest'], st['peaks'],
                                     st['peak_ref'])
    b_off, win, tsys, out_dir = st['b_off'], st['win'], st['tsys'], st['out']
    multi, tag = st['multi'], st['tag']
    inst = st['inst'] or dict(NO_INSTRUMENT, tsys=tsys)

    pairs = []
    for s in sess:
        pairs += make_pairs(s, b_off)

    n = ref['freq'].size
    freq = ref['freq']
    df = freq[1] - freq[0]
    v0 = C * (rest - freq) / rest          # radio velocity before LSR corr.
    dv = C * df / rest                     # km/s per channel (signed)
    good = channel_mask(n, inst)
    if a.inject and not peaks:
        sys.exit('ERROR: --inject needs peak velocities for the target in '
                 'known_targets.json.')

    vref = float(np.mean([p['vsys'] for p in pairs]))
    v = v0 + vref

    for p in pairs:
        ta = tsys * (p['on_data'] - p['off_data']) / p['off_data']
        vk = v0 + p['vsys']                # this pair's own LSR axis
        if a.inject:
            sig = INJ_FWHM / 2.3548
            for vp, _ in peaks:
                ta = ta + a.inject / 1e3 * np.exp(-0.5 * ((vk - vp) / sig) ** 2)
        p['ta'] = np.where(good, ta, np.nan)
        d = p['vsys'] - vref
        p['shift'] = int(np.rint(d / dv))
        p['resid'] = d - p['shift'] * dv
        p['w'] = 1.0 / (1.0 / p['t_on'] + 1.0 / p['t_off'])
        p['rms_th'] = tsys * np.sqrt((1.0 / p['t_on'] + 1.0 / p['t_off'])
                                     / abs(df))
        p['rms'] = baseline(vk, p['ta'], win)[4]

    print(f'\nChannel width: {abs(df) / 1e3:.2f} kHz = {abs(dv):.2f} km/s; '
          f'reference LSR correction (mean of the ON): {vref:+.2f} km/s')
    print(f'\n  {"session":32s} {"ON":>15s} {"OFF":>15s} {"VLSRcor":>8s} '
          f'{"shift":>5s} {"resid":>6s} {"rms":>5s} {"radiom":>6s} {"ratio":>5s}')
    for p in pairs:
        print(f'  {p["session"]:32s} {stamp(p["on_name"]):>15s} '
              f'{stamp(p["off_name"]):>15s} {p["vsys"]:+8.2f} '
              f'{p["shift"]:+5d} {p["resid"]:+6.2f} {p["rms"] * 1e3:5.0f} '
              f'{p["rms_th"] * 1e3:6.0f} {p["rms"] / p["rms_th"]:5.2f}')
    print('  (shift: whole channels applied; resid: remaining misalignment, '
          'km/s; rms and radiometer in mK, each pair on its own)')

    ta_s, wsum = stack(pairs, n)
    res, base, line, spikes, rms = baseline(v, ta_s, win)
    rms_th = tsys / np.sqrt(abs(df) * sum(p['w'] for p in pairs))
    t_on = sum(p['t_on'] for p in pairs)
    t_off = sum(p['t_off'] for p in pairs)
    fin = np.isfinite(res)
    i_pk = np.argmax(np.where(fin & line, res, -np.inf))
    labels = [s['label'] for s in sess]
    dates = [s['utc'][0][:10].replace('-', '') for s in sess]

    print(f'\nSUM: {len(pairs)} pairs from {len(sess)} session(s); '
          f'ON {t_on:.0f} s, OFF {t_off:.0f} s')
    for i in np.flatnonzero(spikes):
        print(f'  ref. channel {i:3d}: {freq[i] / 1e6:.4f} MHz, {v[i]:7.1f} '
              f'km/s, {res[i] * 1e3:+6.0f} mK  (> 5 sigma outside the line '
              'window, excluded from the noise)')
    print(f'Noise outside the line window: {rms * 1e3:.0f} mK '
          f'(radiometer, Tsys {tsys:.0f} K: {rms_th * 1e3:.0f} mK, '
          f'ratio {rms / rms_th:.2f})')
    print(f'Highest channel in {win[0]:g}..{win[1]:g} km/s: '
          f'{res[i_pk] * 1e3:.0f} mK at {v[i_pk]:.1f} km/s '
          f'({res[i_pk] / rms:.1f} sigma)')
    if a.inject:
        print(f'\n*** TEST: synthetic peaks of {a.inject:.0f} mK injected ***')
        for vp, lab in peaks:
            near = np.flatnonzero(fin & (np.abs(v - vp) <= 2.0))
            i = near[np.argmax(res[near])]
            print(f'  {lab} at {vp:+.2f} km/s: recovered {res[i] * 1e3:.0f} '
                  f'mK at {v[i]:+.2f} km/s ({res[i] / rms:.1f} sigma)')

    # ---------------- figure ----------------
    if not multi:
        groups = [([p], stamp(p['on_name'])) for p in pairs]
        gkind = 'each ON/OFF pair'
    else:
        groups = []
        for s, d in zip(labels, dates):
            gp = [p for p in pairs if p['session'] == s]
            groups.append((gp, f'{d} ({len(gp)} pairs)'))
        gkind = 'each session'

    s_ch = v[1] - v[0]                     # km/s per reference channel

    def v2ch(x):
        return (np.asarray(x) - v[0]) / s_ch

    def ch2v(c):
        return v[0] + s_ch * np.asarray(c)

    def v2f(x):
        return (freq[0] + df * v2ch(x)) / 1e6

    def f2v(fm):
        return ch2v((np.asarray(fm) * 1e6 - freq[0]) / df)

    def ref_span(c_lo, c_hi):
        js = [c - p['shift'] for p in pairs for c in (c_lo, c_hi)]
        return sorted((float(ch2v(min(js) - 0.5)), float(ch2v(max(js) + 0.5))))

    vlo, vhi = v[fin].min(), v[fin].max()
    dc_lo, dc_hi = n // 2 - inst['dc_half'], n // 2 + inst['dc_half']
    dc_v = ref_span(dc_lo, dc_hi) if inst['mask_dc'] else None
    spur_v = ref_span(*inst['spur_watch']) if inst['spur_watch'] else None
    pad = 10 if a.channel_axis else 6

    def decorate(ax):
        ax.axvspan(*win, color='lightyellow', zorder=0)
        if dc_v:
            ax.axvspan(*dc_v, color='0.85', zorder=0)
        for vp, _ in peaks:
            ax.axvline(vp, color='tab:green', ls='--', lw=1)
        ax.set_xlim(vlo, vhi)
        if a.channel_axis:
            if spur_v:
                ax.axvspan(*spur_v, color='#ffc8c8', zorder=0)
            sec = ax.secondary_xaxis('top', functions=(v2ch, ch2v))
            sec.xaxis.set_major_locator(MultipleLocator(32))
            sec.xaxis.set_minor_locator(MultipleLocator(8))
            sec.set_xlabel('Channel (reference grid, 0-based)', fontsize=11)

    def channel_ruler(ax):
        """Every channel number along the bottom of a panel (SVG only)."""
        trb = mtransforms.blended_transform_factory(ax.transData,
                                                    ax.transAxes)
        js = np.flatnonzero((v >= vlo) & (v <= vhi))
        ax.vlines(v[js], 0.0, 0.010, transform=trb, color='0.5', lw=0.2)
        for j in js:
            ax.text(v[j], 0.012, str(j), transform=trb, rotation=90,
                    fontsize=RULER_FONT, color='0.3', ha='center',
                    va='bottom')

    fig, (a1, a2) = plt.subplots(2, 1, figsize=FIG_SIZE)

    gres = []
    for gp, lab in groups:
        s_, _ = stack(gp, n)
        r_, _, _, _, rms_g = baseline(v, s_, win)
        gres.append((r_, lab, rms_g))
    step = 6.0 * np.median([g[2] for g in gres]) * 1e3
    tr1 = mtransforms.blended_transform_factory(a1.transAxes, a1.transData)
    decorate(a1)
    for k, (r_, lab, _) in enumerate(gres):
        y0 = -k * step
        a1.axhline(y0, color='0.8', lw=0.6)
        a1.plot(v, r_ * 1e3 + y0, color='black', lw=0.8, drawstyle='steps-mid')
        a1.text(1.012, y0, lab, transform=tr1, va='center', fontsize=11)
    a1.set_yticks([])
    a1.set_ylabel(f'T_A, offset {step:.0f} mK per row')
    a1.set_title(f'{gkind.capitalize()} (baseline removed), on the common '
                 'LSR axis', pad=pad)

    rp = res * 1e3
    decorate(a2)
    a2.plot(v, rp, color='black', lw=1, drawstyle='steps-mid')
    if spikes.any():
        a2.plot(v[spikes], rp[spikes], 'o', color='tab:red', ms=5,
                label='> 5 sigma outside the line window (interference?)')
        for i in np.flatnonzero(spikes):
            a2.annotate(f'ch {i}', (v[i], rp[i]), xytext=(6, 0),
                        textcoords='offset points', color='tab:red',
                        fontsize=10, va='center', ha='left')
        a2.legend(loc='lower right')
    a2.axhline(0, color='0.6', lw=0.8)
    for sgn in (1, -1):
        a2.axhline(sgn * rms * 1e3, color='0.6', lw=0.8, ls=':')
    ytop = np.nanmax(rp)
    for k, (vp, lab) in enumerate(peaks):
        a2.text(vp, ytop * (0.95 - 0.1 * k), ' ' + lab, color='tab:green',
                fontsize=11, va='top')
    tr2 = mtransforms.blended_transform_factory(a2.transData, a2.transAxes)
    if dc_v:
        a2.text(dc_v[1], 0.06, f' band centre, ch {dc_lo}-{dc_hi} (masked)',
                transform=tr2, color='0.35', fontsize=10, va='bottom',
                ha='left')
    if a.channel_axis and spur_v:
        a2.text(spur_v[0], 0.06, f'ch {inst["spur_watch"][0]}-'
                f'{inst["spur_watch"][1]}: candidate spur (watched, not '
                'masked) ', transform=tr2, color='darkred', fontsize=10,
                va='bottom', ha='right')
    a2.set_xlabel(f'LSR velocity for {rest / 1e6:.3f} MHz (km/s)')
    a2.set_ylabel('T_A (mK, nominal)')
    a2.set_title(f'Sum: rms {rms * 1e3:.0f} mK (radiometer '
                 f'{rms_th * 1e3:.0f} mK, ratio {rms / rms_th:.2f}); '
                 'yellow: expected line range, excluded from the baseline',
                 pad=pad)

    fax = a2.secondary_xaxis(FREQ_AXIS_POS, functions=(v2f, f2v))
    fax.xaxis.set_major_locator(MultipleLocator(0.1))
    fax.xaxis.set_minor_locator(MultipleLocator(0.02))
    fax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    fax.tick_params(axis='x', which='both', colors=FREQ_COLOR)
    for sp in fax.spines.values():
        sp.set_color(FREQ_COLOR)
    if multi:
        fax.set_xlabel(f'Sky frequency as recorded (MHz), reference grid: '
                       f'LSR correction {vref:+.2f} km/s', color=FREQ_COLOR)
    else:
        fax.set_xlabel('Sky frequency as recorded (MHz)', color=FREQ_COLOR)

    span = dates[0] if not multi else f'{dates[0]} - {dates[-1]}'
    title = (f'{facility} {instrument} - {target or "unnamed target"}, rest '
             f'{rest / 1e6:.3f} MHz, position switching: {len(pairs)} ON/OFF '
             f'pairs, {len(sess)} session(s), {span}')
    if a.inject:
        title = (f'TEST - SYNTHETIC PEAKS OF {a.inject:.0f} mK INJECTED - '
                 + title)
    fig.suptitle(title, fontsize=15, color='red' if a.inject else 'black')

    notes = [
        f'T_A = Tsys (ON - OFF) / OFF for each pair, nominal Tsys = '
        f'{tsys:.0f} K ({inst["note"]}); LSR correction of each ON from the '
        f'.fro (Doppler/vlsr_correction_mps); pairs aligned by whole-channel '
        f'shifts (no interpolation, residual <= {abs(dv) / 2:.2f} km/s) and '
        f'summed with radiometer weights; baseline: polynomial of degree '
        f'{BASE_DEG}.',
        (f'Dashed: peak velocities from {peak_ref}. ' if peaks else '')
        + (f'Grey: band centre, channels {dc_lo}-{dc_hi} (masked). ' if dc_v else '')
        + 'Red dots: channels above 5 sigma outside the line window, labelled '
          'with their channel (0-based, reference grid). SVG version: every '
          'channel number along the bottom of each panel (zoom in).']
    if a.channel_axis and spur_v:
        notes.append(f'Pink: channels {inst["spur_watch"][0]}-'
                     f'{inst["spur_watch"][1]}, candidate internal spur '
                     '(highlighted only, not masked).')
    # Notes as wide as the plots: wrap them to the width of the lower
    # panel, then lay out again with room for the resulting lines.
    note_lines = notes
    for _ in range(3):
        notes_h = 0.05 + len(note_lines) * NOTE_FONT * 1.45 / 72  # inches
        fig.tight_layout(rect=(0, notes_h / fig.get_figheight(), 0.94, 0.96))
        box = a2.get_position()
        width_pt = box.width * fig.get_figwidth() * 72
        chars = max(60, int(width_pt / (NOTE_FONT * NOTE_CHAR_EM)))
        new_lines = [ln for nt in notes for ln in textwrap.wrap(nt, chars)]
        if len(new_lines) == len(note_lines):
            note_lines = new_lines
            break
        note_lines = new_lines
    box = a2.get_position()
    fig.text(box.x0 + box.width / 2, 0.006, '\n'.join(note_lines),
             ha='center', va='bottom', fontsize=NOTE_FONT)

    os.makedirs(out_dir, exist_ok=True)
    if multi:
        stem = (f'{facility}_{tag}_stack_{dates[0]}-{dates[-1]}_'
                f'{len(sess)}sess')
    else:
        stem = f'{ref["label"]}_onoff'
    if a.inject:
        stem += '_INJECTED_TEST'
    out = os.path.join(out_dir, stem)
    fig.savefig(out + '.png', dpi=120)
    channel_ruler(a1)                      # SVG only: added after the PNG
    channel_ruler(a2)
    fig.savefig(out + '.svg')
    np.savez(out + '.npz', v_lsr_kms=v, ta_k=res, ta_stack_k=ta_s,
             baseline_k=base, weight_sum=wsum, rms_k=rms,
             rms_radiometer_k=rms_th, tsys_k=tsys, vref_kms=vref,
             rest_frequency_hz=rest, line_window_kms=np.array(win),
             off_glat_deg=b_off, injected_mk=a.inject,
             sessions=np.array(labels),
             fro_files=np.array([os.path.basename(s['path']) for s in sess]),
             on_spectra=np.array([p['on_name'] for p in pairs]),
             off_spectra=np.array([p['off_name'] for p in pairs]))
    for ext in ('.png', '.svg', '.npz'):
        print('Saved:', out + ext)
    if not a.no_open:
        viewer = next((v for v in SVG_VIEWERS if shutil.which(v)), 'xdg-open')
        try:
            subprocess.Popen([viewer, out + '.svg'],
                             stdout=subprocess.DEVNULL,
                             stderr=subprocess.DEVNULL,
                             start_new_session=True)
            print(f'Opened with {viewer}:', out + '.svg')
        except FileNotFoundError:
            print('No viewer available: open the SVG above by hand.')


if __name__ == '__main__':
    main()
