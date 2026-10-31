#!/usr/bin/env python3
"""
plot_salsa_oh_onoff.py - SALSA position-switched spectra at the OH 1612 MHz
line, for one session or for the sum of several sessions (days).

For every ON/OFF pair:  T_A = Tsys * (ON - OFF) / OFF, on an LSR velocity
axis for the OH rest frequency, using the VELOSYS of that pair's ON spectrum
only. The pairs are aligned on a common LSR axis by whole-channel shifts (no
interpolation, so the noise statistics are untouched; residual misalignment
<= half a channel, ~0.45 km/s) and summed with radiometer weights. A
polynomial baseline is fitted outside the expected line window and the
measured noise is compared with the radiometer equation.

Usage:
  plot_salsa_oh_onoff.py [options] <session folder> [<session folder> ...]

Options:
  --off-lat DEG    galactic latitude of the OFF position (default 12)
  --out DIR        output folder (default: the session folder for a single
                   session; <parent>/SALSA_OH127_1612_stack for several)
  --channel-axis   add a channel axis on top of both panels and highlight
                   the watched spur channels
  --inject MK      TEST ONLY: add two synthetic peaks of MK mK at the Wolak
                   et al. 2013 velocities to every pair before stacking, to
                   check alignment and sensitivity. Never use for results.

Outputs: PNG, SVG and NPZ. The SVG also carries every channel number along
the bottom of each panel, in tiny print, readable when zoomed in.

The sum panel has a second axis (blue) with the sky frequency as recorded by
SALSA (CRVAL1/CDELT1). With several sessions it refers to the common
reference grid, i.e. to a pair whose VELOSYS equals the reference value.

ON = files whose GLAT differs from the OFF latitude by more than 1 deg.
Inside each session the files are taken in name (= time) order; with equal
numbers of ON and OFF the k-th ON is paired with the k-th OFF, otherwise each
ON with the nearest OFF in time (a warning is printed: the radiometer
estimate then assumes independent OFFs and is optimistic).

Channels: 0-based. The "reference grid" is the channel grid of the sum; a
channel c of a pair lands on reference channel c - shift.
"""
import argparse
import glob
import os
import sys
import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.ticker import FormatStrFormatter, MultipleLocator
from astropy.io import fits

warnings.filterwarnings('ignore')

C = 299792.458              # km/s
F0_OH = 1612.231e6          # Hz, OH ground state F=1-2 satellite line
TSYS = 300.0                # K, SALSA nominal (approximate, see SALSA docs)
EDGE = 40                   # channels masked at each band edge
DC_HALF = 1                 # channels masked each side of the band centre
SPUR_CHANNELS = ()          # fixed internal spurs to MASK (0-based channel
                            # numbers); empty until confirmed on a second
                            # session
SPUR_WATCH = (320, 322)     # channel range only HIGHLIGHTED with
                            # --channel-axis (candidate spur, 2026-10-28)
LINE_WIN = (-75.0, -35.0)   # km/s LSR, excluded from the baseline fit
PEAKS = ((-65.9, 'blue peak'), (-43.65, 'red peak'))  # Wolak et al. 2013
INJ_FWHM = 1.5              # km/s, width of the synthetic test peaks
BASE_DEG = 3                # baseline polynomial degree
RULER_FONT = 2.5            # pt, channel numbers in the SVG ruler
FREQ_AXIS_POS = -0.22       # frequency axis position, axes fraction
FREQ_COLOR = 'tab:blue'     # colour of the frequency axis
STACK_DIR = 'SALSA_OH127_1612_stack'

plt.rcParams['font.size'] = 12


def load(folder):
    rows = []
    for p in sorted(glob.glob(os.path.join(folder, '*.fits'))):
        with fits.open(p) as f:
            h = f[0].header
            d = np.asarray(f[0].data, dtype=float).ravel()
        rows.append(dict(
            name=os.path.basename(p), data=d,
            glat=float(h['CRVAL3']),
            f0=float(h['CRVAL1']), df=float(h['CDELT1']),
            pix=float(h['CRPIX1']), t=float(h['OBSTIME']),
            vsys=float(h.get('VELOSYS', 0.0)) / 1000.0,
            el=float(h.get('ELEVATIO', np.nan))))
    return rows


def make_pairs(folder, b_off):
    """Read one session folder and return (rows, ON/OFF pairs)."""
    rows = load(folder)
    if not rows:
        sys.exit(f'ERROR: no FITS files in {folder}')
    session = os.path.basename(os.path.normpath(folder))
    print(f'\nSession {session}')
    print(f'  {"file":36s} {"GLAT":>6s} {"El":>6s} {"t (s)":>6s} '
          f'{"VELOSYS":>8s}  role')
    for k, r in enumerate(rows):
        r['k'] = k
        r['off'] = abs(r['glat'] - b_off) <= 1.0
        print(f'  {r["name"]:36s} {r["glat"]:6.1f} {r["el"]:6.1f} '
              f'{r["t"]:6.0f} {r["vsys"]:8.2f}  {"OFF" if r["off"] else "ON"}')
    on = [r for r in rows if not r['off']]
    off = [r for r in rows if r['off']]
    if not on or not off:
        sys.exit(f'ERROR: {session}: need both ON and OFF '
                 f'({len(on)} ON, {len(off)} OFF)')
    if len(on) == len(off):
        match = list(zip(on, off))
    else:
        print(f'  WARNING: {len(on)} ON and {len(off)} OFF: each ON paired '
              'with the nearest OFF in time (radiometer estimate optimistic)')
        match = [(a, min(off, key=lambda b: abs(b['k'] - a['k'])))
                 for a in on]
    return rows, [dict(session=session, on=a, off=b) for a, b in match]


def channel_mask(n):
    good = np.ones(n, bool)
    good[:EDGE] = good[-EDGE:] = False
    dc = n // 2                                 # band centre (LO leakage)
    good[dc - DC_HALF:dc + DC_HALF + 1] = False
    for c in SPUR_CHANNELS:
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


def baseline(v, ta):
    """Polynomial baseline outside the line window; robust noise.

    Narrow interference is searched for only outside the line window, so a
    real (or injected) line is never flagged as interference."""
    n = ta.size
    valid = np.isfinite(ta)
    line = (v > LINE_WIN[0]) & (v < LINE_WIN[1])
    fit = valid & ~line
    x = (np.arange(n) - n / 2) / (n / 2)
    base = np.polyval(np.polyfit(x[fit], ta[fit], BASE_DEG), x)
    res = ta - base
    mad = 1.4826 * np.median(np.abs(res[fit] - np.median(res[fit])))
    spikes = fit & (np.abs(np.nan_to_num(res)) > 5.0 * mad)
    rms = np.std(res[fit & ~spikes])
    return res, base, line, spikes, rms


def stamp(name):
    return name.split('-')[-1].replace('.fits', '')


def main():
    ap = argparse.ArgumentParser(
        description='SALSA OH 1612 MHz ON/OFF spectra, one or more sessions.')
    ap.add_argument('folders', nargs='+', help='session folder(s) with FITS')
    ap.add_argument('--off-lat', type=float, default=12.0,
                    help='GLAT of the OFF position (default 12)')
    ap.add_argument('--out', help='output folder')
    ap.add_argument('--channel-axis', action='store_true',
                    help='channel axis on top and highlighted watched spur '
                         'channels')
    ap.add_argument('--inject', type=float, default=0.0,
                    help='TEST ONLY: synthetic peak amplitude in mK')
    a = ap.parse_args()

    folders = [os.path.abspath(f) for f in a.folders]
    for f in folders:
        if not os.path.isdir(f):
            sys.exit(f'ERROR: not a folder: {f}')

    allrows, pairs = [], []
    for f in folders:
        rows, pp = make_pairs(f, a.off_lat)
        allrows += rows
        pairs += pp

    ref = allrows[0]
    n = ref['data'].size
    same = all(r['data'].size == n and r['f0'] == ref['f0']
               and r['df'] == ref['df'] and r['pix'] == ref['pix']
               for r in allrows)
    if not same:
        sys.exit('ERROR: frequency axes differ between files '
                 '(all sessions must use the same settings)')

    df = ref['df']
    freq = ref['f0'] + (np.arange(n) + 1 - ref['pix']) * df
    v0 = C * (F0_OH - freq) / F0_OH        # radio velocity before LSR corr.
    dv = C * df / F0_OH                    # km/s per channel (signed)
    good = channel_mask(n)

    # Reference LSR correction: mean VELOSYS of the ON spectra (with a single
    # session this reproduces the axis of the first version exactly).
    vref = float(np.mean([p['on']['vsys'] for p in pairs]))
    v = v0 + vref

    for p in pairs:
        on, off = p['on']['data'], p['off']['data']
        ta = TSYS * (on - off) / off
        vk = v0 + p['on']['vsys']          # this pair's own LSR axis
        if a.inject:
            sig = INJ_FWHM / 2.3548
            for vp, _ in PEAKS:
                ta = ta + a.inject / 1e3 * np.exp(-0.5 * ((vk - vp) / sig) ** 2)
        p['ta'] = np.where(good, ta, np.nan)
        d = p['on']['vsys'] - vref
        p['shift'] = int(np.rint(d / dv))
        p['resid'] = d - p['shift'] * dv
        p['w'] = 1.0 / (1.0 / p['on']['t'] + 1.0 / p['off']['t'])
        p['rms_th'] = TSYS * np.sqrt((1.0 / p['on']['t']
                                      + 1.0 / p['off']['t']) / df)
        p['rms'] = baseline(vk, p['ta'])[4]

    print(f'\nChannel width: {df / 1e3:.2f} kHz = {abs(dv):.2f} km/s; '
          f'reference LSR correction (mean VELOSYS of the ON): {vref:+.2f} km/s')
    print(f'\n  {"session":28s} {"ON":>15s} {"OFF":>15s} {"VELOSYS":>8s} '
          f'{"shift":>5s} {"resid":>6s} {"rms":>5s} {"radiom":>6s} {"ratio":>5s}')
    for p in pairs:
        print(f'  {p["session"]:28s} {stamp(p["on"]["name"]):>15s} '
              f'{stamp(p["off"]["name"]):>15s} '
              f'{p["on"]["vsys"]:+8.2f} {p["shift"]:+5d} {p["resid"]:+6.2f} '
              f'{p["rms"] * 1e3:5.0f} {p["rms_th"] * 1e3:6.0f} '
              f'{p["rms"] / p["rms_th"]:5.2f}')
    print('  (shift: whole channels applied; resid: remaining misalignment, '
          'km/s; rms and radiometer in mK, each pair on its own)')

    ta_s, wsum = stack(pairs, n)
    res, base, line, spikes, rms = baseline(v, ta_s)
    rms_th = TSYS / np.sqrt(df * sum(p['w'] for p in pairs))
    t_on = sum(p['on']['t'] for p in pairs)
    t_off = sum(p['off']['t'] for p in pairs)
    fin = np.isfinite(res)
    sel = fin & line
    i_pk = np.argmax(np.where(sel, res, -np.inf))

    sessions = []
    for p in pairs:
        if p['session'] not in sessions:
            sessions.append(p['session'])

    print(f'\nSUM: {len(pairs)} pairs from {len(sessions)} session(s); '
          f'ON {t_on:.0f} s, OFF {t_off:.0f} s')
    for i in np.flatnonzero(spikes):
        print(f'  ref. channel {i:3d}: {freq[i] / 1e6:.4f} MHz, {v[i]:7.1f} '
              f'km/s, {res[i] * 1e3:+6.0f} mK  (> 5 sigma outside the line '
              'window, excluded from the noise)')
    print(f'Noise outside the line window: {rms * 1e3:.0f} mK '
          f'(radiometer, Tsys {TSYS:.0f} K: {rms_th * 1e3:.0f} mK, '
          f'ratio {rms / rms_th:.2f})')
    print(f'Highest channel in {LINE_WIN[0]:.0f}..{LINE_WIN[1]:.0f} km/s: '
          f'{res[i_pk] * 1e3:.0f} mK at {v[i_pk]:.1f} km/s '
          f'({res[i_pk] / rms:.1f} sigma)')
    if a.inject:
        print(f'\n*** TEST: synthetic peaks of {a.inject:.0f} mK injected ***')
        for vp, lab in PEAKS:
            near = np.flatnonzero(fin & (np.abs(v - vp) <= 2.0))
            i = near[np.argmax(res[near])]
            print(f'  {lab} at {vp:+.2f} km/s: recovered {res[i] * 1e3:.0f} '
                  f'mK at {v[i]:+.2f} km/s ({res[i] / rms:.1f} sigma)')

    # ---------------- figure ----------------
    if len(sessions) == 1:
        groups = [([p], stamp(p['on']['name'])) for p in pairs]
        gkind = 'each ON/OFF pair'
    else:
        groups = []
        for s in sessions:
            gp = [p for p in pairs if p['session'] == s]
            groups.append((gp, f'{s.split("_")[-1]} ({len(gp)} pairs)'))
        gkind = 'each session'

    # Reference channel <-> LSR velocity <-> sky frequency (all linear).
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
        """Velocity span covered by pair channels c_lo..c_hi in all pairs."""
        js = [c - p['shift'] for p in pairs for c in (c_lo, c_hi)]
        return sorted((float(ch2v(min(js) - 0.5)), float(ch2v(max(js) + 0.5))))

    vlo, vhi = v[fin].min(), v[fin].max()
    spur_v = ref_span(*SPUR_WATCH)
    dc_lo, dc_hi = n // 2 - DC_HALF, n // 2 + DC_HALF
    dc_v = ref_span(dc_lo, dc_hi)
    pad = 10 if a.channel_axis else 6

    def decorate(ax):
        ax.axvspan(*LINE_WIN, color='lightyellow', zorder=0)
        ax.axvspan(*dc_v, color='0.85', zorder=0)
        for vp, _ in PEAKS:
            ax.axvline(vp, color='tab:green', ls='--', lw=1)
        ax.set_xlim(vlo, vhi)
        if a.channel_axis:
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

    fig, (a1, a2) = plt.subplots(2, 1, figsize=(21, 10))

    gres = []
    for gp, lab in groups:
        s_, _ = stack(gp, n)
        r_, _, _, _, rms_g = baseline(v, s_)
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
    for s in (1, -1):
        a2.axhline(s * rms * 1e3, color='0.6', lw=0.8, ls=':')
    ytop = np.nanmax(rp)
    for k, (vp, lab) in enumerate(PEAKS):
        a2.text(vp, ytop * (0.95 - 0.1 * k), ' ' + lab, color='tab:green',
                fontsize=11, va='top')
    tr2 = mtransforms.blended_transform_factory(a2.transData, a2.transAxes)
    a2.text(dc_v[1], 0.06, f' band centre, ch {dc_lo}-{dc_hi} (masked)',
            transform=tr2, color='0.35', fontsize=10, va='bottom', ha='left')
    if a.channel_axis:
        a2.text(spur_v[0], 0.06, f'ch {SPUR_WATCH[0]}-{SPUR_WATCH[1]}: '
                'candidate spur (watched, not masked) ', transform=tr2,
                color='darkred', fontsize=10, va='bottom', ha='right')
    a2.set_xlabel('LSR velocity for OH 1612.231 MHz (km/s)')
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
    if len(sessions) == 1:
        fax.set_xlabel('Sky frequency as recorded by SALSA (MHz)',
                       color=FREQ_COLOR)
    else:
        fax.set_xlabel(f'Sky frequency as recorded by SALSA (MHz), reference '
                       f'grid: VELOSYS {vref:+.2f} km/s', color=FREQ_COLOR)

    span = (sessions[0].split('_')[-1] if len(sessions) == 1 else
            f'{sessions[0].split("_")[-1]} - {sessions[-1].split("_")[-1]}')
    title = (f'SALSA - OH127.8+0.0, OH 1612 MHz, position switching: '
             f'{len(pairs)} ON/OFF pairs, {len(sessions)} session(s), {span}')
    if a.inject:
        title = (f'TEST - SYNTHETIC PEAKS OF {a.inject:.0f} mK INJECTED - '
                 + title)
    fig.suptitle(title, fontsize=15, color='red' if a.inject else 'black')

    notes = [
        f'T_A = Tsys (ON - OFF) / OFF for each pair, nominal Tsys = '
        f'{TSYS:.0f} K; LSR correction: VELOSYS of each ON; pairs aligned by '
        f'whole-channel shifts (no interpolation, residual <= '
        f'{abs(dv) / 2:.2f} km/s) and summed with radiometer weights; '
        f'baseline: polynomial of degree {BASE_DEG}.',
        f'Dashed: peak velocities from Wolak et al. 2013. Grey: band centre, '
        f'channels {dc_lo}-{dc_hi} (masked). Red dots: channels above 5 sigma '
        'outside the line window, labelled with their channel (0-based, '
        'reference grid). SVG version: every channel number along the bottom '
        'of each panel (zoom in).']
    if a.channel_axis:
        notes.append(f'Pink: channels {SPUR_WATCH[0]}-{SPUR_WATCH[1]}, '
                     'candidate internal spur seen on 2026-10-28 '
                     '(highlighted only, not masked).')
    fig.text(0.5, 0.006, '\n'.join(notes), ha='center', va='bottom',
             fontsize=10)
    fig.tight_layout(rect=(0, 0.005 + 0.02 * len(notes), 0.94, 0.96))

    if a.out:
        out_dir = os.path.abspath(a.out)
    elif len(folders) == 1:
        out_dir = folders[0]
    else:
        out_dir = os.path.join(os.path.dirname(folders[0]), STACK_DIR)
    os.makedirs(out_dir, exist_ok=True)
    if len(folders) == 1:
        stem = 'salsa_oh_onoff'
    else:
        stem = (f'salsa_oh_stack_{sessions[0].split("_")[-1]}-'
                f'{sessions[-1].split("_")[-1]}_{len(sessions)}sess')
    if a.inject:
        stem += '_INJECTED_TEST'
    out = os.path.join(out_dir, stem)
    fig.savefig(out + '.png', dpi=120)
    channel_ruler(a1)                      # SVG only: added after the PNG
    channel_ruler(a2)
    fig.savefig(out + '.svg')
    np.savez(out + '.npz', v_lsr_kms=v, ta_k=res, ta_stack_k=ta_s,
             baseline_k=base, weight_sum=wsum, rms_k=rms,
             rms_radiometer_k=rms_th, tsys_k=TSYS, vref_kms=vref,
             injected_mk=a.inject, sessions=np.array(sessions),
             on_files=np.array([p['on']['name'] for p in pairs]),
             off_files=np.array([p['off']['name'] for p in pairs]))
    for ext in ('.png', '.svg', '.npz'):
        print('Saved:', out + ext)


if __name__ == '__main__':
    main()
