#!/bin/bash
# Launcher for gbt_sdfits_to_fro.py: graphical FITS file selection
# via zenity, then conversion. Terminal stays open at the end so the
# user can read the output.

FITS_FILE=$(zenity --file-selection \
    --title="Select a GBT SDFITS file" \
    --filename="$HOME/FRO/GBT/" \
    --file-filter="FITS files | *.fits *.FITS" \
    --file-filter="All files | *")

if [ -z "$FITS_FILE" ]; then
    echo "No file selected. Exiting."
    read -p "Press Enter to close..."
    exit 0
fi

echo "Converting: $FITS_FILE"
echo "----------------------------------------"
python3 "$HOME/FRO/Parsers/GBT/gbt_sdfits_to_fro.py" "$FITS_FILE"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
