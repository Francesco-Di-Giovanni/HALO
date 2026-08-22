# effelsberg_to_fro.py
# Parser: MBFITS (Effelsberg 100m, EDD backend) -> FRO HDF5 format v17
# One .fro file per baseband, all subscans merged.
# DATA has NUSEFEED polarization channels per integration (FEBEPAR:
# USEFEED/BESECTS show these are backend sections of a single physical
# feed, not separate beams) -> written as NUSEFEED spectra per integration.
# DATAPAR PHASE (signal/reference switching phase, not an ON/CAL/OFF state)
# is stored per spectrum in flag_cal; separation is left to analysis.
# Reference observation: Holmberg1, scan 6723, FEBE P217mm-EDD.
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy

import os
import sys
import glob
import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from datetime import datetime, timezone

sys.path.insert(0, os.path.expanduser('~/FRO/FRO_System'))
from fro_format_v17 import create_fro_v17_structure, append_raw_spectrum, read_fro_summary

# Observatory: Effelsberg 100m (same values as ebhis_to_fro.py)
OBS_NAME = 'Effelsberg'
OBS_LAT = 50.5248    # deg N
OBS_LON = 6.8836     # deg E
OBS_ALT = 369.0      # m

# Position-switching classification: LONGOFF/LATOFF (deg) below this
# radius from the nominal source position are treated as ON-source.
ON_OFF_OFFSET_THRESHOLD_DEG = 1.0

_base = os.path.expanduser('~/FRO/Parsers/Effelsberg/Surveys/Holmberg1')
MBFITS_DIR = os.path.join(_base, 'Holmberg1_MBFITS_Profiles',
                           'EFFBG_2026-06-23_6723_02-26_P217mm-EDD')
OUTPUT_DIR = os.path.join(_base, 'Holmberg1_FRO')


def get_febe_name(mbfits_dir):
    scan = fits.open(os.path.join(mbfits_dir, 'SCAN.fits'))
    febe = scan[1].data['FEBE'][0]
    scan.close()
    return febe


def get_top_level_metadata(mbfits_dir):
    """Reads the primary header of GROUPING.fits (present in every MBFITS
    obs directory) for OBJECT/TELESCOP/SCANTYPE/SCANMODE/PROJID.
    Falls back to generic defaults if the file is missing."""
    top_path = os.path.join(mbfits_dir, 'GROUPING.fits')
    meta = {'object': 'unknown', 'telescop': OBS_NAME,
            'scantype': 'unknown', 'scanmode': 'unknown', 'projid': ''}
    if os.path.isfile(top_path):
        hdu = fits.open(top_path)
        hdr = hdu[0].header
        meta['object']   = str(hdr.get('OBJECT', meta['object'])).strip()
        meta['telescop'] = str(hdr.get('TELESCOP', meta['telescop'])).strip()
        meta['scantype'] = str(hdr.get('SCANTYPE', meta['scantype'])).strip()
        meta['scanmode'] = str(hdr.get('SCANMODE', meta['scanmode'])).strip()
        meta['projid']   = str(hdr.get('PROJID', meta['projid'])).strip()
        hdu.close()
    return meta


def get_subscan_dirs(mbfits_dir):
    dirs = sorted([
        d for d in glob.glob(os.path.join(mbfits_dir, '*'))
        if os.path.isdir(d) and os.path.basename(d).isdigit()
    ], key=lambda d: int(os.path.basename(d)))
    return dirs


def get_basebands(subscan_dir, febe):
    pattern = os.path.join(subscan_dir, febe + '-ARRAYDATA-*.fits')
    files = sorted(glob.glob(pattern))
    basebands = [int(os.path.basename(f).split('-')[-1].replace('.fits', ''))
                 for f in files]
    return basebands


def get_frequency_axis(arraydata_header, n_channels):
    crpx = arraydata_header['1CRPX2F']
    crvl = arraydata_header['1CRVL2F']
    cd   = arraydata_header['11CD2F']
    channels = np.arange(1, n_channels + 1)
    freq_hz = crvl + (channels - crpx) * cd
    return freq_hz


def mjd_to_datetime(mjd):
    t = Time(mjd, format='mjd', scale='utc')
    return t.to_datetime(timezone=timezone.utc)


def compute_galactic(ra_deg, dec_deg):
    coord = SkyCoord(ra=ra_deg*u.deg, dec=dec_deg*u.deg, frame='icrs')
    return coord.galactic.l.deg, coord.galactic.b.deg


def parse_mbfits(mbfits_dir, output_dir):
    obs_name = os.path.basename(mbfits_dir.rstrip('/'))
    febe = get_febe_name(mbfits_dir)
    meta = get_top_level_metadata(mbfits_dir)
    print('FEBE:', febe)
    print('Object:', meta['object'], '| Scan type:', meta['scantype'], '/', meta['scanmode'])

    subscan_dirs = get_subscan_dirs(mbfits_dir)
    print('Subscans found:', len(subscan_dirs))

    basebands = get_basebands(subscan_dirs[0], febe)
    print('Basebands:', basebands)

    os.makedirs(output_dir, exist_ok=True)

    for bb in basebands:
        print('\n── Baseband', bb, '──────────────────────────────')

        ad_path = os.path.join(subscan_dirs[0], febe + '-ARRAYDATA-' + str(bb) + '.fits')
        ad_hdu = fits.open(ad_path)
        ad_hdr = ad_hdu[1].header
        n_channels = ad_hdr['CHANNELS']
        n_feeds = int(ad_hdr.get('NUSEFEED', 1))
        freq_hz = get_frequency_axis(ad_hdr, n_channels)
        freq_mhz = freq_hz / 1e6
        rest_freq_hz = float(ad_hdr.get('RESTFREQ', 0.0))
        transition = str(ad_hdr.get('TRANSITI', 'unknown')).strip()
        bandwidth_hz = float(ad_hdr.get('BANDWID', 0.0))
        date_obs = ad_hdr.get('DATE-OBS', 'unknown')
        scan_num = int(ad_hdr.get('SCANNUM', 0))
        ad_hdu.close()

        print('  Channels:', n_channels)
        print('  Feeds (polarization channels):', n_feeds)
        print('  Freq range:', round(freq_mhz[0], 3), '–', round(freq_mhz[-1], 3), 'MHz')
        print('  Transition:', transition)
        print('  Rest freq:', round(rest_freq_hz/1e9, 6), 'GHz')

        out_name = obs_name + '_bb' + str(bb) + '.fro'
        out_path = os.path.join(output_dir, out_name)

        pol_labels = ['pol' + str(i + 1) for i in range(n_feeds)]

        notes = (
            'MBFITS format (Effelsberg 100m, EDD backend). Observation: ' + obs_name + '. '
            'Baseband ' + str(bb) + ' of ' + str(len(basebands)) + '. '
            'Object: ' + meta['object'] + '. Scan type/mode: ' + meta['scantype'] + '/' + meta['scanmode'] + '. '
            'Transition: ' + transition + ' at ' + str(round(rest_freq_hz/1e9, 6)) + ' GHz. '
            'Date-obs: ' + str(date_obs) + '. '
            'DATA array has ' + str(n_feeds) + ' feed(s)/polarization channel(s) per integration '
            '(FEBEPAR USEFEED/BESECTS show these are backend polarization sections of a single '
            'physical feed, not separate beams; POLTY is undefined ("N") in the source FEBEPAR, so '
            'channels are labeled generically ' + ', '.join(pol_labels) + '). '
            'flag_cal stores the raw MBFITS DATAPAR PHASE value (1 or 2) for each spectrum: an '
            'internal signal/reference switching phase (RA/DEC unchanged between phases, WOBUSED=F '
            'in the source data), not an ON/CAL/OFF calibration state. Separation of signal vs '
            'reference phases is left to the analysis stage. '
            'flag_track marks position-switched ON/OFF from DATAPAR LONGOFF/LATOFF: 1 when the '
            'offset from the nominal source position is below ' + str(ON_OFF_OFFSET_THRESHOLD_DEG) + ' deg (on source), '
            '0 otherwise (reference sky position; this dataset alternates OFF-ON-OFF-ON by subscan). '
            'Frequency axis from ARRAYDATA header (LSR rest frame). '
            'Site: Effelsberg 100m (lat=' + str(OBS_LAT) + ', lon=' + str(OBS_LON) + ', alt=' + str(OBS_ALT) + 'm).'
        )

        create_fro_v17_structure(
            out_path,
            frequency_axis_hz   = freq_hz,
            facility            = meta['telescop'] if meta['telescop'] else OBS_NAME,
            instrument          = febe,
            observation_id      = str(scan_num),
            survey              = 'Holmberg1',
            dataset             = obs_name,
            provider            = 'MPIfR Effelsberg',
            import_date         = datetime.now(timezone.utc).isoformat(),
            observatory_name    = OBS_NAME,
            latitude_deg        = OBS_LAT,
            longitude_deg       = OBS_LON,
            elevation_m         = OBS_ALT,
            receiver            = febe,
            center_frequency_hz = rest_freq_hz,
            bandwidth_hz        = bandwidth_hz,
            pointing_source     = 'measured',
            target              = meta['object'],
            scan_id             = scan_num,
            subscan_type        = meta['scantype'] + '/' + meta['scanmode'],
            signal              = 'unknown',
            schedule_name       = meta['projid'],
            cal_polarizations   = pol_labels,
            cal_mark_temp_k     = [float('nan')] * n_feeds,
            attenuation_db      = [float('nan')] * n_feeds,
            notes               = notes,
        )
        print('  Created:', out_path)

        total_spectra = 0

        for ss_dir in subscan_dirs:
            ss_num = int(os.path.basename(ss_dir))

            dp_path = os.path.join(ss_dir, febe + '-DATAPAR.fits')
            dp_hdu = fits.open(dp_path)
            dp = dp_hdu[1].data
            n_integ = len(dp)
            dp_hdu.close()

            ad_path = os.path.join(ss_dir, febe + '-ARRAYDATA-' + str(bb) + '.fits')
            ad_hdu = fits.open(ad_path)
            spectra = ad_hdu[1].data['DATA']  # shape (n_integ, n_feeds, n_channels)
            ad_hdu.close()

            print('  Subscan', ss_num, ':', n_integ, 'integrations, spectra shape', spectra.shape)

            for i in range(n_integ):
                mjd      = float(dp['MJD'][i])
                integtim = float(dp['INTEGTIM'][i])
                ra_deg   = float(dp['RA'][i])
                dec_deg  = float(dp['DEC'][i])
                az_deg   = float(dp['AZIMUTH'][i])
                el_deg   = float(dp['ELEVATIO'][i])
                phase    = float(dp['PHASE'][i])
                longoff  = float(dp['LONGOFF'][i])
                latoff   = float(dp['LATOFF'][i])
                glon, glat = compute_galactic(ra_deg, dec_deg)
                utc_dt   = mjd_to_datetime(mjd)

                # Position-switched ON/OFF: LONGOFF/LATOFF near zero means
                # pointing on the nominal source position (ON), a large
                # offset means a reference sky position (OFF). This flips
                # by subscan for Holmberg1 (subscans 1,3 = OFF; 2,4 = ON).
                on_source = (longoff**2 + latoff**2) ** 0.5 < ON_OFF_OFFSET_THRESHOLD_DEG
                flag_track = 1 if on_source else 0

                for feed_idx in range(n_feeds):
                    spec = spectra[i, feed_idx].astype(np.float32)
                    append_raw_spectrum(
                        out_path,
                        raw_spectrum       = spec,
                        utc_time           = utc_dt,
                        integration_time_s = integtim,
                        az_deg             = az_deg,
                        elev_deg           = el_deg,
                        ra_deg             = ra_deg,
                        dec_deg            = dec_deg,
                        glon_deg           = glon,
                        glat_deg           = glat,
                        polarization       = pol_labels[feed_idx],
                        flag_track         = flag_track,
                        flag_cal           = phase,
                    )
                    total_spectra += 1

        print('  Total spectra written:', total_spectra)
        read_fro_summary(out_path)


if __name__ == '__main__':
    if len(sys.argv) > 2:
        parse_mbfits(sys.argv[1], sys.argv[2])
    elif len(sys.argv) > 1:
        parse_mbfits(sys.argv[1], OUTPUT_DIR)
    else:
        parse_mbfits(MBFITS_DIR, OUTPUT_DIR)
