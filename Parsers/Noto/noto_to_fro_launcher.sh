#!/bin/bash
# noto_to_fro_launcher.sh
FITS_FILE=$(zenity --file-selection     --title="Select Noto FITS file"     --filename="$HOME/FRO/Parsers/Noto/Surveys/"     --file-filter="FITS files (*.fits *.fit) | *.fits *.fit"     --file-filter="All files | *"     2>/dev/null)
if [ -z "$FITS_FILE" ]; then
    echo "No file selected. Exiting."
    exit 0
fi
echo "Converting: $FITS_FILE"
echo "----------------------------------------"
env -u LD_LIBRARY_PATH python3 "$HOME/FRO/Parsers/Noto/noto_to_fro.py" "$FITS_FILE"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
