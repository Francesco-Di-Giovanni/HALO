#!/bin/bash
# wsrt_to_fro_launcher.sh
# Supports both MeasurementSet (.ms directory) and FITS input

MODE=$(zenity --list \
    --title="WSRT / Apertif → HALO" \
    --text="Tipo di input:" \
    --column="Tipo" \
    "FITS cube / UVFITS (.fits)" \
    "MeasurementSet (.ms directory)" \
    2>/dev/null)

if [ -z "$MODE" ]; then
    exit 0
fi

if [ "$MODE" = "MeasurementSet (.ms directory)" ]; then
    INPUT=$(zenity --file-selection --directory \
        --title="Seleziona directory MeasurementSet (.ms)" \
        --filename="$HOME/FRO/Parsers/WSRT/Surveys/" \
        2>/dev/null)
else
    INPUT=$(zenity --file-selection \
        --title="Seleziona FITS cube o UVFITS" \
        --filename="$HOME/FRO/Parsers/WSRT/Surveys/" \
        --file-filter="FITS files | *.fits *.fit" \
        --file-filter="All files | *" \
        2>/dev/null)
fi

if [ -z "$INPUT" ]; then
    echo "Nessun file selezionato. Uscita."
    exit 0
fi

echo "Converting: $INPUT"
echo "----------------------------------------"
env -u LD_LIBRARY_PATH python3 "$HOME/FRO/Parsers/WSRT/wsrt_to_fro.py" "$INPUT"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
