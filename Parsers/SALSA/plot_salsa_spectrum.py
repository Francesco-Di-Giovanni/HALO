#!/usr/bin/env python3
import sys, os, warnings
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

warnings.filterwarnings('ignore')
C_KMS = 299792.458
HI_MHZ = 1420.405752

def load_spectrum(path):
    hdul = fits.open(path)
    for hdu in hdul:
        if hdu.data is not None and hdu.data.squeeze().ndim == 1:
            data = hdu.data.squeeze().astype(float)
            h = hdu.header
            if 'CRVAL1' in h and 'CDELT1' in h:
                n = len(data)
                crpix = h.get('CRPIX1', 1)
                freqs = (h['CRVAL1'] + (np.arange(n) - crpix + 1) * h['CDELT1']) / 1e6
                vlsr = C_KMS * (HI_MHZ - freqs) / HI_MHZ
                return vlsr, data, os.path.basename(path)
    return None, None, None

if len(sys.argv) < 2:
    print("Usage: plot_salsa_spectrum.py <file.fits>")
    sys.exit(1)

vlsr, data, title = load_spectrum(sys.argv[1])
if vlsr is None:
    print("Cannot read spectrum")
    sys.exit(1)

plt.figure(figsize=(10, 4))
plt.plot(vlsr, data, 'b-', linewidth=0.8)
plt.xlabel('VLSR (km/s)')
plt.ylabel('Amplitude')
plt.title(title)
plt.xlim(-300, 300)
plt.axvline(0, color='r', linewidth=0.5, linestyle='--')
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.show()
