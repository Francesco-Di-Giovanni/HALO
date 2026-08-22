#!/bin/bash
# lab_to_fro_launcher.sh
# Desktop launcher wrapper for lab_to_fro.py

cd "$HOME/FRO/Parsers/LAB"
echo "Starting LAB to FRO converter..."
echo "----------------------------------------"
env -u LD_LIBRARY_PATH python3 "$HOME/FRO/Parsers/LAB/lab_to_fro.py"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
