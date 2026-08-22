#!/usr/bin/env python3
# rename_salsa_fits.py
#
# Renames SALSA FITS files adding GLon/GLat prefix computed from
# the pointing coordinates in the FITS header.
# Format: G{lon:06.2f}{lat:+.2f}_{original_name}.fits
# Example: G140.00+05.00_SALSA-torre-20260710T090103.fits
#
# Handles two SALSA FITS formats:
#   - New format: BinTable HDU with AZIMUTH/ALTITUDE keywords
#   - Old format: PRIMARY only with RA/DEC WCS (CRVAL2/CRVAL3)
#
# Usage:
#   python3 rename_salsa_fits.py <folder>           # dry-run
#   python3 rename_salsa_fits.py <folder> --apply   # actually rename

import os
import sys
import argparse
import warnings

from astropy.io import fits
from astropy.coordinates import SkyCoord, AltAz, EarthLocation
from astropy.time import Time
import astropy.units as u

# Onsala Observatory (SALSA Torre)
ONSALA = EarthLocation(lat=57.3953*u.deg, lon=11.9255*u.deg, height=10*u.m)


def get_glon_glat(fits_path):
    """Extract GLon/GLat from a SALSA FITS file.
    Handles both formats:
    - New format: BinTable HDU with AZIMUTH/ALTITUDE keywords
    - Old format: PRIMARY only with RA/DEC WCS keywords (CRVAL2/CRVAL3)
    """
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        hdul = fits.open(fits_path)

    # Try new format: any HDU with AZIMUTH/ALTITUDE keyword
    hdr = None
    for hdu in hdul:
        if "AZIMUTH" in hdu.header or "AZ" in hdu.header:
            hdr = hdu.header
            break

    if hdr is not None:
        # New format: Az/El -> AltAz -> Galactic
        az    = float(hdr.get("AZIMUTH",  hdr.get("AZ",  0.0)))
        el    = float(hdr.get("ALTITUDE", hdr.get("EL",  0.0)))
        t_str = hdr.get("DATE-OBS", hdr.get("DATE", None))
        hdul.close()
        if t_str is None:
            return None, None
        t = Time(t_str, format="isot", scale="utc")
        altaz_frame = AltAz(obstime=t, location=ONSALA)
        altaz = SkyCoord(az=az*u.deg, alt=el*u.deg, frame=altaz_frame)
        gal = altaz.icrs.galactic
        return float(gal.l.deg), float(gal.b.deg)

    # Old format: RA/DEC WCS in PRIMARY header
    hdr = hdul[0].header
    hdul.close()
    ra  = hdr.get("CRVAL2", None)
    dec = hdr.get("CRVAL3", None)
    if ra is None or dec is None:
        return None, None
    gal = SkyCoord(ra=float(ra)*u.deg,
                   dec=float(dec)*u.deg, frame="icrs").galactic
    return float(gal.l.deg), float(gal.b.deg)


def make_prefix(glon, glat):
    """Build the GLon/GLat prefix string.
    Rounds to nearest 5 degrees for survey grid alignment.
    Example: G140+05 or G095-10
    """
    glon_r = int(round(glon / 5.0)) * 5
    glat_r = int(round(glat / 5.0)) * 5
    return f"G{glon_r:03d}{glat_r:+03d}"


def rename_salsa_folder(folder, apply=False):
    """Rename all SALSA FITS files in folder with GLon/GLat prefix."""
    files = sorted([f for f in os.listdir(folder)
                    if f.endswith(".fits") and "SALSA" in f
                    and not f.startswith("G")])  # skip already renamed

    print(f"Found {len(files)} files to rename (apply={apply})\n")
    errors = []
    renamed = 0
    for fname in files:
        old_path = os.path.join(folder, fname)
        try:
            glon, glat = get_glon_glat(old_path)
            if glon is None:
                print(f"  SKIP (no coords): {fname}")
                continue
            prefix = make_prefix(glon, glat)
            new_name = f"{prefix}_{fname}"
            new_path = os.path.join(folder, new_name)
            print(f"  {fname}")
            print(f"    -> {new_name}")
            if apply:
                os.rename(old_path, new_path)
            renamed += 1
        except Exception as e:
            print(f"  ERROR {fname}: {e}")
            errors.append(fname)

    print(f"\nDone. {renamed} files {'renamed' if apply else 'would be renamed'}, "
          f"{len(errors)} errors.")
    if errors:
        print("Files with errors:")
        for e in errors:
            print(f"  {e}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Rename SALSA FITS files with GLon/GLat prefix.")
    parser.add_argument("folder", help="Folder with SALSA FITS files")
    parser.add_argument("--apply", action="store_true",
                        help="Actually rename files (default: dry-run)")
    args = parser.parse_args()
    rename_salsa_folder(args.folder, apply=args.apply)