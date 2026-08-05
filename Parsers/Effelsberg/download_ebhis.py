# download_ebhis.py
# Download HI profiles from the Bonn HI Survey server (EBHIS/HI4PI)
# for a grid of galactic coordinates.
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy

import os
import time
import urllib.request
import urllib.parse
import re

# ── Grid settings — asked interactively at runtime ─────────────────────────────
GLON_MIN = float(input('GLon min (deg): '))
GLON_MAX = float(input('GLon max (deg): '))
GLAT_MIN = float(input('GLat min (deg): '))
GLAT_MAX = float(input('GLat max (deg): '))
STEP     = float(input('Step (deg) [2.5]: ') or '2.5')
BEAM     = 0.200  # degrees, matches EBHIS native resolution

# Auto-build directory names from grid coordinates
GRID_TAG   = f'G{GLON_MIN:03.0f}-{GLON_MAX:03.0f}_B{GLAT_MIN:+03.0f}{GLAT_MAX:+03.0f}'
SURVEY_DIR = os.path.expanduser(f'~/FRO/Parsers/Effelsberg/Surveys/EBHIS_Survey_{GRID_TAG}')
OUTPUT_DIR = os.path.join(SURVEY_DIR, f'EBHIS_Profiles_{GRID_TAG}')
# ───────────────────────────────────────────────────────────────────────────────
BASE_URL = 'https://www.astro.uni-bonn.de/hisurvey/AllSky_profiles/'
FORM_URL = BASE_URL + 'index.php'
DOWNLOAD_URL = BASE_URL + 'download.php'
DELAY = 2.0   # seconds between requests, be polite to the server


def download_profile(glon, glat, output_dir):
    """Download one HI profile from the Bonn HI survey server."""

    filename = 'EBHIS_G%06.2f%+06.2f.txt' % (glon, glat)
    filepath = os.path.join(output_dir, filename)

    if os.path.exists(filepath):
        print('  SKIP (exists):', filename)
        return True

    # Step 1: POST to generate the data
    form_data = urllib.parse.urlencode({
        'coordinates': 'lb',
        'ral': '%.2f' % glon,
        'decb': '%.2f' % glat,
        'beam': '%.3f' % BEAM,
        'ebhis': 'ebhis',
        'search': 'Search',
    }).encode()

    try:
        req = urllib.request.Request(FORM_URL, data=form_data, method='POST')
        resp = urllib.request.urlopen(req, timeout=30)
        html = resp.read().decode('utf-8', errors='replace')

        # Step 2: find download link
        match = re.search(r'href="(download\.php\?[^"]+)"', html)
        if not match:
            print('  ERROR: no download link found for', filename)
            return False

        dl_url = BASE_URL + match.group(1)

        # Step 3: download the data
        urllib.request.urlretrieve(dl_url, filepath)
        size = os.path.getsize(filepath)
        print('  OK:', filename, '(%d bytes)' % size)
        return True

    except Exception as e:
        print('  ERROR:', filename, str(e))
        return False


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    glons = []
    g = GLON_MIN
    while g <= GLON_MAX + 0.01:
        glons.append(round(g, 1))
        g += STEP

    glats = []
    g = GLAT_MIN
    while g <= GLAT_MAX + 0.01:
        glats.append(round(g, 1))
        g += STEP

    total = len(glons) * len(glats)
    print('Download grid: GLon %.1f-%.1f, GLat %.1f-%.1f, step %.1f deg' %
          (GLON_MIN, GLON_MAX, GLAT_MIN, GLAT_MAX, STEP))
    print('Total pointings:', total)
    print('Output:', OUTPUT_DIR)
    print()

    done = 0
    errors = 0

    for glon in glons:
        for glat in glats:
            done += 1
            print('[%d/%d] GLon=%.1f GLat=%.1f' % (done, total, glon, glat))
            ok = download_profile(glon, glat, OUTPUT_DIR)
            if not ok:
                errors += 1
            time.sleep(DELAY)

    print()
    print('Done: %d downloaded, %d errors' % (done - errors, errors))


if __name__ == '__main__':
    main()