# mbfits_to_fro.py
# Parser: MBFITS (APEX/Effelsberg) → FRO HDF5 format v17
# One .fro file per baseband, all subscans merged.
# Reference file: APEX-48698-2011-08-04-T-087.F-0001-2011/ (format v1.63)
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

# ── User settings — asked interactively at runtime ─────────────────────────────
_base = os.path.expanduser('~/FRO/Parsers/APEX/')
OBS_TAG  = input('Observation tag (e.g. 12CO_345GHz): ').strip()
SURVEY_DIR = os.path.join(_base, f'MBFITS_Obs_{OBS_TAG}')
MBFITS_DIR = os.path.join(SURVEY_DIR, f'MBFITS_Profiles_{OBS_TAG}')
OUTPUT_DIR = os.path.join(SURVEY_DIR, f'MBFITS_FRO_{OBS_TAG}')
# ───────────────────────────────────────────────────────────────────────────────

def get_febe_name(mbfits_dir):
    scan = fits.open(os.path.join(mbfits_dir, 'SCAN.fits'))
    febe = scan[1].data['FEBE'][0]
    scan.close()
    return febe


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
    obs_name = os.path.basename(mbfits_dir)
    febe = get_febe_name(mbfits_dir)
    print('FEBE:', febe)

    subscan_dirs = get_subscan_dirs(mbfits_dir)
    print('Subscans found:', len(subscan_dirs))

    basebands = get_basebands(subscan_dirs[0], febe)
    print('Basebands:', basebands)

    for bb in basebands:
        print('\n── Baseband', bb, '──────────────────────────────')

        ad_path = os.path.join(subscan_dirs[0], febe + '-ARRAYDATA-' + str(bb) + '.fits')
        ad_hdu = fits.open(ad_path)
        ad_hdr = ad_hdu[1].header
        n_channels = ad_hdr['CHANNELS']
        freq_hz = get_frequency_axis(ad_hdr, n_channels)
        freq_mhz = freq_hz / 1e6
        rest_freq_hz = float(ad_hdr.get('RESTFREQ', 0.0))
        transition = ad_hdr.get('TRANSITI', 'unknown').strip()
        bandwidth_hz = float(ad_hdr.get('BANDWID', 0.0))
        date_obs = ad_hdr.get('DATE-OBS', 'unknown')
        scan_num = int(ad_hdr.get('SCANNUM', 0))
        ad_hdu.close()

        print('  Channels:', n_channels)
        print('  Freq range:', round(freq_mhz[0], 3), '–', round(freq_mhz[-1], 3), 'MHz')
        print('  Transition:', transition)
        print('  Rest freq:', round(rest_freq_hz/1e9, 6), 'GHz')

        out_name = obs_name + '_bb' + str(bb) + '.fro'
        out_path = os.path.join(output_dir, out_name)

        notes = ('MBFITS format v1.63. Observation: ' + obs_name + '. '
                 'Baseband ' + str(bb) + ' of ' + str(len(basebands)) + '. '
                 'Transition: ' + transition + ' at ' + str(round(rest_freq_hz/1e9, 6)) + ' GHz. '
                 'Date-obs: ' + date_obs + '. '
                 'Frequency axis from ARRAYDATA header (LSR rest frame). '
                 'Site: APEX (lat=-22.9534, lon=-67.7592, alt=5105m).')

        create_fro_v17_structure(
            out_path,
            frequency_axis_hz = freq_hz,
            facility          = 'APEX' if 'APEX' in obs_name else 'Effelsberg',
            instrument        = febe,
            observation_id    = str(scan_num),
            survey            = 'MBFITS sample',
            dataset           = obs_name,
            provider          = 'NASA FITS Registry',
            import_date       = datetime.now(timezone.utc).isoformat(),
            bandwidth_hz      = bandwidth_hz,
            pointing_source   = 'measured',
            target            = transition,
            notes             = notes,
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
            spectra = ad_hdu[1].data['DATA']
            ad_hdu.close()

            print('  Subscan', ss_num, ':', n_integ, 'integrations, spectra shape', spectra.shape)

            for i in range(n_integ):
                mjd       = float(dp['MJD'][i])
                integtim  = float(dp['INTEGTIM'][i])
                ra_deg    = float(dp['RA'][i])
                dec_deg   = float(dp['DEC'][i])
                az_deg    = float(dp['AZIMUTH'][i])
                el_deg    = float(dp['ELEVATIO'][i])
                glon, glat = compute_galactic(ra_deg, dec_deg)
                utc_dt    = mjd_to_datetime(mjd)
                spec      = spectra[i].flatten().astype(np.float32)

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
                    polarization       = 'unknown',
                    flag_track         = -1,
                    flag_cal           = 0,
                )
                total_spectra += 1

        print('  Total spectra written:', total_spectra)
        read_fro_summary(out_path)


if __name__ == '__main__':
    import sys
    if len(sys.argv) > 1:
        parse_mbfits(sys.argv[1], OUTPUT_DIR)
    else:
        parse_mbfits(MBFITS_DIR, OUTPUT_DIR)