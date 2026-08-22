#!/bin/bash
# Launcher for APEX MBFITS plots — static + interactive
# FRO/HALO Project — Francesco Di Giovanni
PLOTS_DIR=~/FRO/Parsers/APEX/Surveys/MBFITS_Obs_12CO_345GHz/MBFITS_Plots_12CO_345GHz

zenity --info \
  --title="APEX MBFITS — 12CO(3-2)" \
  --text="Dataset: APEX-48698-2011-08-04\nStrumento: HET345-XFFTS2\nGLon: 208.75 - 208.99 deg\nGLat: -19.83 - -19.39 deg\nSpettri: 81 per banda base (bb1+bb2)" \
  --width=400 2>/dev/null &

# Apri solo i PNG statici
for f in "$PLOTS_DIR"/*.png; do
    [ -f "$f" ] && xdg-open "$f" &
done

# Apri il viewer interattivo (spettri con cursore)
cd ~/FRO/Parsers/APEX
env -u LD_LIBRARY_PATH python3 ~/FRO/Parsers/APEX/view_apex_interactive.py
