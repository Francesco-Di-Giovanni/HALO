#!/bin/bash
# parkes_gass_to_fro_launcher.sh
# Desktop launcher wrapper for parkes_gass_to_fro.py:
# graphical FITS file selection via zenity

FITS_FILE=$(zenity --file-selection \
    --title="Select Parkes GASS SDFITS file" \
    --filename="$HOME/FRO/Parsers/Parkes/Surveys/" \
    --file-filter="FITS files (*.fits *.fit) | *.fits *.fit" \
    --file-filter="All files | *" \
    2>/dev/null)

if [ -z "$FITS_FILE" ]; then
    echo "No file selected. Exiting."
    exit 0
fi

echo "Converting: $FITS_FILE"
echo "----------------------------------------"
env -u LD_LIBRARY_PATH python3 "$HOME/FRO/Parsers/Parkes/parkes_gass_to_fro.py" "$FITS_FILE"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
