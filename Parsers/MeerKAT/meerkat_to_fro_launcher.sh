#!/bin/bash
# meerkat_to_fro_launcher.sh
# Supports both MeasurementSet (.ms directory) and FITS cube input

MODE=$(zenity --list \
    --title="MeerKAT → FRO" \
    --text="Tipo di input:" \
    --column="Tipo" \
    "FITS cube (.fits)" \
    "MeasurementSet (.ms directory)" \
    2>/dev/null)

if [ -z "$MODE" ]; then
    exit 0
fi

if [ "$MODE" = "MeasurementSet (.ms directory)" ]; then
    INPUT=$(zenity --file-selection --directory \
        --title="Seleziona directory MeasurementSet (.ms)" \
        --filename="$HOME/FRO/Parsers/MeerKAT/Surveys/" \
        2>/dev/null)
else
    INPUT=$(zenity --file-selection \
        --title="Seleziona FITS cube" \
        --filename="$HOME/FRO/Parsers/MeerKAT/Surveys/" \
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
env -u LD_LIBRARY_PATH python3 "$HOME/FRO/Parsers/MeerKAT/meerkat_to_fro.py" "$INPUT"
echo "----------------------------------------"
echo "Done."
read -p "Press Enter to close..."
