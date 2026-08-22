#!/bin/bash
FITS_FILE=$(zenity --file-selection --title="Select FAST FITS file" --file-filter="FITS files | *.fits *.fit" 2>/dev/null)
if [ -z "$FITS_FILE" ]; then
    echo "No file selected. Exiting."
    exit 1
fi
KY_TABLE=$(zenity --file-selection --title="Select KY RA/Dec table (Table-mjd-radec.txt) or Cancel to skip" --file-filter="Text files | *.txt" 2>/dev/null)
OUT_DIR=$(dirname "$FITS_FILE")
cd /home/franz/FRO/FAST
if [ -z "$KY_TABLE" ]; then
    echo "No KY table selected — coordinates will be NaN."
    python3 /home/franz/FRO/Parsers/FAST/fast_fits_to_fro.py "$FITS_FILE" "$OUT_DIR"
else
    python3 /home/franz/FRO/Parsers/FAST/fast_fits_to_fro.py "$FITS_FILE" "$OUT_DIR" "$KY_TABLE"
fi
echo "Press Enter to close..."
read
