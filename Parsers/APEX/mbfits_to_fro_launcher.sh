#!/bin/bash
# Launcher for mbfits_to_fro.py: graphical directory selection
# via zenity, then conversion. Terminal stays open at the end.
MBFITS_DIR=$(zenity --file-selection     --title="Select MBFITS directory"     --filename="$HOME/FRO/Parsers/Effelsberg/"     --directory)
if [ -z "$MBFITS_DIR" ]; then
    echo "No directory selected. Exiting."
    read -p "Press Enter to close..."
    exit 0
fi
echo "Converting: $MBFITS_DIR"
echo "----------------------------------------"
env -u LD_LIBRARY_PATH python3 /home/franz/FRO/Parsers/Effelsberg/mbfits_to_fro.py "$MBFITS_DIR"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
