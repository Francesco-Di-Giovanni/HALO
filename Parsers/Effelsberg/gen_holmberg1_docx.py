import copy
from docx import Document

SRT_IT = "/home/franz/FRO/Parsers/SRT/SRT_Data_Structure_IT.docx"
SRT_EN = "/home/franz/FRO/Parsers/SRT/SRT_Data_Structure_EN.docx"
OUT_IT = "/home/franz/FRO/Parsers/Effelsberg/Surveys/Holmberg1/Holmberg1_Data_Structure_IT.docx"
OUT_EN = "/home/franz/FRO/Parsers/Effelsberg/Surveys/Holmberg1/Holmberg1_Data_Structure_EN.docx"

# ---- contenuto: (tipo, testo) ----
# T=title, S=sub, H=heading, P=paragrafo, B=paragrafo tutto grassetto
IT = [
("T","Struttura Dati Effelsberg/Holmberg1 e Modalità Operativa"),
("S","FRO / HALO (Hydrogen Atomic Line Observatory) — Francesco Di Giovanni, Bolzano, Italia"),
("S","Versione documento: 1.0 — 9 agosto 2026"),
("H","1. Panoramica dell'Osservazione"),
("P","Il parser effelsberg_to_fro.py legge osservazioni in formato MBFITS del radiotelescopio Effelsberg da 100 m (MPIfR, Bonn) e le converte nel formato FRO HDF5 (Hierarchical Data Format 5) v17. È lo stesso formato dati MBFITS già gestito dal parser APEX (mbfits_to_fro.py), con cui condivide la logica di lettura dell'asse di frequenza; le differenze specifiche di Effelsberg sono descritte alla sezione 7. Il set di dati di riferimento per validazione è lo scan 6723 (progetto 02-26) sulla galassia nana Holmberg I, acquisito con il ricevitore/backend P217mm-EDD (banda 21 cm)."),
("H","2. Telescopio e Strumento"),
("P","Tipo: radiotelescopio single-dish da 100 m a montatura altazimutale. Sito: Effelsberg, Germania — quota 369 m s.l.m. Coordinate: lat 50.5248°N, lon 6.8836°E (stessi valori già usati dal parser EBHIS in questa stessa cartella). Operatore: MPIfR (Max-Planck-Institut für Radioastronomie), Bonn. FEBE (combinazione frontend-backend) usata in questo scan: P217mm-EDD, ricevitore a 21 cm accoppiato al backend digitale EDD (Effelsberg Direct Digitization). Formato dati: MBFITS v1.67 (stesso standard MBFITS del parser APEX)."),
("H","3. Formato MBFITS e Differenze rispetto ad APEX"),
("P","La struttura directory è identica a quella APEX: un file SCAN.fits a livello di osservazione (colonna FEBE), una sottocartella numerata per ogni subscan contenente FEBE-ARRAYDATA-N.fits (i dati) e FEBE-DATAPAR.fits (puntamento e stato per ogni integrazione), più un FEBE-FEBEPAR.fits e un GROUPING.fits a livello di osservazione. L'asse di frequenza si ricava dalle stesse parole chiave header (1CRPX2F, 1CRVL2F, 11CD2F, CHANNELS) usate da mbfits_to_fro.py."),
("P","Due differenze strutturali reali rispetto al set di riferimento APEX (NUSEFEED=1, un solo canale per integrazione):"),
("B","NUSEFEED=2 (canali di polarizzazione)."),
("P","In questo scan ogni riga di ARRAYDATA contiene 2 canali (forma dei dati: integrazioni x 2 x canali), non un fascio spaziale separato: FEBEPAR mostra USEFEED=[1,1] e BESECTS=[1,2], cioè le due sezioni del backend applicate allo stesso feed fisico. POLTY nel FEBEPAR sorgente è indefinito (\"N\"), quindi il parser etichetta i due canali genericamente pol1/pol2 invece di H/V o L/R. Ogni integrazione produce quindi 2 spettri .fro."),
("B","Colonna PHASE in DATAPAR."),
("P","DATAPAR contiene una colonna PHASE (valori 1/2) assente nell'uso che ne fa il parser APEX. Nei dati di Holmberg1 questa alterna riga per riga a RA/DEC costanti entro un subscan (WOBUSED=F), quindi non è uno spostamento in cielo ma una fase di switching segnale/riferimento interna al backend. Il parser scrive il valore PHASE grezzo nel campo Quality/flag_cal di ogni spettro; la separazione segnale/riferimento è lasciata alla fase di analisi."),
("H","4. Dataset di Riferimento: Holmberg1 (scan 6723)"),
("B","Perché questo target:"),
("P","Holmberg I (UGC 5139) è una galassia nana irregolare del gruppo di M81, nota in letteratura per il suo contenuto di idrogeno neutro (riga HI a 21 cm). Lo scan 6723 (progetto 02-26, osservatore/operatore \"SAINTONGE/WS\" come registrato nel file di log dell'osservazione) è un puntamento ON/OFF sulla sorgente."),
("B","Cosa è stato verificato:"),
("P","Il parser è stato eseguito end-to-end sui dati grezzi reali. Le 4 sottoscansioni contengono rispettivamente 980, 898, 984 e 896 integrazioni (3758 totali), ciascuna scritta come 2 spettri (pol1/pol2) per un totale di 7516 spettri nel file .fro prodotto. Il confronto tra i valori LONGOFF/LATOFF dei 4 subscan conferma un pattern di position-switching OFF-ON-OFF-ON (subscan 1 e 3 con offset ~3.09° dalla posizione nominale = OFF/riferimento; subscan 2 e 4 con offset ~0° = ON/sorgente), coerente con il flag_track scritto dal parser (si veda sezione 7). Non è stata ancora eseguita alcuna analisi spettrale (baseline, ricerca riga) né generato alcun plot: questo documento descrive la sola validazione dell'importazione dati grezzi."),
("P","Caratteristiche tecniche: 262144 canali per spettro, larghezza di banda 400 MHz, risoluzione in frequenza 1525.879 Hz/canale, intervallo di frequenza 1200.002-1600.0 MHz. RESTFREQ nell'header è 1.400000 GHz esatto: è un valore nominale di sintonia, non la frequenza di riposo precisa della riga HI (1420.405752 MHz); la conversione in velocità LSR andrà fatta con la frequenza di riposo reale della riga, non con RESTFREQ."),
("H","5. Provenienza dei Dati"),
("P","A differenza di APEX/COHRS/EBHIS, questo non è un archivio pubblico scaricabile: è un'osservazione propria (progetto Effelsberg 02-26), archiviata localmente in ~/FRO/Parsers/Effelsberg/Surveys/Holmberg1/Holmberg1 Raw/ (archivio Holmberg1_Ef_data.tag.gz, ~3.4 GB decompresso, più le immagini Holmberg1.png/Holmberg1_zoom.png e observation.log del turno osservativo). L'archivio MBFITS estratto risiede in Surveys/Holmberg1/Holmberg1_MBFITS_Profiles/."),
("H","6. Modalità Operativa"),
("B","Passo 1 — Estrarre l'archivio MBFITS."),
("P","Decomprimere Holmberg1_Ef_data.tag.gz nella cartella Holmberg1_MBFITS_Profiles/ (l'archivio è in formato tar.gz standard nonostante l'estensione .tag.gz)."),
("B","Passo 2 — Eseguire il parser."),
("P","env -u LD_LIBRARY_PATH python3 ~/FRO/Parsers/Effelsberg/effelsberg_to_fro.py [MBFITS_DIR] [OUTPUT_DIR]. Senza argomenti usa i percorsi di default per Holmberg1. Un file .fro viene scritto per baseband in Holmberg1_FRO/."),
("B","Passo 3 — Verificare e archiviare."),
("P","Controllare il riepilogo stampato da read_fro_summary() (conteggio spettri, range RA/DEC/GLON/GLAT, righe non-tracking). Eseguire l'alias fro-backup per rispecchiare ~/FRO/ in ~/FRO_Backup/."),
("H","7. Parser: effelsberg_to_fro.py"),
("P","Posizione: ~/FRO/Parsers/Effelsberg/. Un file .fro per baseband (in questo scan: 1 solo baseband), tutte le sottoscansioni unite nello stesso file. Metadati di osservazione (OBJECT, TELESCOP, SCANTYPE, SCANMODE, PROJID) letti dall'header primario (HDU0) di GROUPING.fits dentro la cartella dell'osservazione — non da un ipotetico file esterno, poiché nell'archivio il file .fits di riepilogo a livello superiore è in realtà un symlink relativo a questo stesso GROUPING.fits. facility='Effelsberg 100m RT' (da TELESCOP), target='HOI' (da OBJECT). polarization='pol1'/'pol2' per i due canali NUSEFEED. flag_cal = valore grezzo PHASE (1.0 o 2.0) di DATAPAR. flag_track = 1 se la distanza (LONGOFF, LATOFF) dalla posizione nominale è inferiore a 1.0° (ON, costante ON_OFF_OFFSET_THRESHOLD_DEG nello script), altrimenti 0 (OFF). Sito, coordinate e note di calibrazione temperatura sono coerenti con ebhis_to_fro.py nella stessa cartella. Niente è hardcoded: percorso MBFITS da argv o default, tutti i parametri di frequenza/polarizzazione/puntamento letti dagli header FITS."),
("H","8. Compatibilità con ezRA"),
("P","La struttura HDF5 interna del formato v17 (/Observatory, /Spectra, /Time, /Pointing) è la stessa attesa da fro_to_ezra.py, quindi i file .fro prodotti da questo parser sono convertibili al formato testuale di ezRA/ezCon.py senza modifiche allo script di conversione."),
("H","9. Riferimenti"),
("P","Specifica formato FRO: fro_format_v17.py (~/FRO/FRO_System/)."),
("P","Formato MBFITS: stesso standard descritto in APEX_MBFITS_Data_Structure_IT_v1.1.docx (~/FRO/Parsers/APEX/), condiviso tra APEX ed Effelsberg."),
("P","MPIfR Effelsberg 100m: Max-Planck-Institut für Radioastronomie, Bonn, Germania."),
("P","Nota: non sono state incluse citazioni bibliografiche specifiche su Holmberg I per evitare di riportare riferimenti non verificati in questa sessione; da integrare se Francesco indica un DOI/paper di riferimento preferito."),
]

EN = [
("T","Effelsberg/Holmberg1 Data Structure and Operating Procedure"),
("S","FRO / HALO (Hydrogen Atomic Line Observatory) — Francesco Di Giovanni, Bolzano, Italy"),
("S","Document version: 1.0 — 9 August 2026"),
("H","1. Observation Overview"),
("P","The parser effelsberg_to_fro.py reads MBFITS-format observations from the Effelsberg 100-m radio telescope (MPIfR, Bonn) and converts them into FRO HDF5 (Hierarchical Data Format 5) v17 format. This is the same MBFITS data format already handled by the APEX parser (mbfits_to_fro.py), with which it shares the frequency-axis reading logic; the Effelsberg-specific differences are described in section 7. The reference dataset for validation is scan 6723 (project 02-26) on the dwarf galaxy Holmberg I, taken with the P217mm-EDD receiver/backend (21 cm band)."),
("H","2. Telescope and Instrument"),
("P","Type: 100-m alt-azimuth single-dish radio telescope. Site: Effelsberg, Germany — altitude 369 m a.s.l. Coordinates: lat 50.5248°N, lon 6.8836°E (same values already used by the EBHIS parser in this same folder). Operator: MPIfR (Max-Planck-Institut für Radioastronomie), Bonn. FEBE (frontend-backend combination) used in this scan: P217mm-EDD, a 21 cm receiver coupled to the EDD (Effelsberg Direct Digitization) digital backend. Data format: MBFITS v1.67 (same MBFITS standard as the APEX parser)."),
("H","3. MBFITS Format and Differences from APEX"),
("P","The directory structure is identical to APEX: a top-level SCAN.fits (FEBE column), one numbered subfolder per subscan containing FEBE-ARRAYDATA-N.fits (data) and FEBE-DATAPAR.fits (per-integration pointing and state), plus an observation-level FEBE-FEBEPAR.fits and GROUPING.fits. The frequency axis is derived from the same header keywords (1CRPX2F, 1CRVL2F, 11CD2F, CHANNELS) used by mbfits_to_fro.py."),
("P","Two real structural differences from the APEX reference set (NUSEFEED=1, a single channel per integration):"),
("B","NUSEFEED=2 (polarization channels)."),
("P","In this scan every ARRAYDATA row holds 2 channels (data shape: integrations x 2 x channels), not a separate spatial beam: FEBEPAR shows USEFEED=[1,1] and BESECTS=[1,2], i.e. the two backend sections applied to the same physical feed. POLTY in the source FEBEPAR is undefined (\"N\"), so the parser labels the two channels generically pol1/pol2 rather than H/V or L/R. Each integration therefore produces 2 .fro spectra."),
("B","PHASE column in DATAPAR."),
("P","DATAPAR contains a PHASE column (values 1/2) not used by the APEX parser. In the Holmberg1 data this alternates row by row at constant RA/DEC within a subscan (WOBUSED=F), so it is not a sky-position switch but an internal signal/reference switching phase of the backend. The parser writes the raw PHASE value into each spectrum's Quality/flag_cal field; separating signal from reference is left to the analysis stage."),
("H","4. Reference Dataset: Holmberg1 (scan 6723)"),
("B","Why this target:"),
("P","Holmberg I (UGC 5139) is a dwarf irregular galaxy in the M81 group, known in the literature for its neutral hydrogen content (21 cm HI line). Scan 6723 (project 02-26, observer/operator \"SAINTONGE/WS\" as recorded in the observation log) is an ON/OFF pointing on the source."),
("B","What was verified:"),
("P","The parser was run end-to-end on the real raw data. The 4 subscans contain 980, 898, 984 and 896 integrations respectively (3758 total), each written as 2 spectra (pol1/pol2) for a total of 7516 spectra in the output .fro file. Comparing LONGOFF/LATOFF across the 4 subscans confirms an OFF-ON-OFF-ON position-switching pattern (subscans 1 and 3 with ~3.09° offset from the nominal position = OFF/reference; subscans 2 and 4 with ~0° offset = ON/source), consistent with the flag_track written by the parser (see section 7). No spectral analysis (baselining, line search) has been performed yet, and no plot has been generated: this document describes raw-data import validation only."),
("P","Technical characteristics: 262144 channels per spectrum, 400 MHz bandwidth, 1525.879 Hz/channel frequency resolution, frequency range 1200.002-1600.0 MHz. The header RESTFREQ is exactly 1.400000 GHz: this is a nominal tuning value, not the precise HI line rest frequency (1420.405752 MHz); LSR velocity conversion should use the true line rest frequency, not RESTFREQ."),
("H","5. Data Provenance"),
("P","Unlike APEX/COHRS/EBHIS, this is not a downloadable public archive: it is a proprietary observation (Effelsberg project 02-26), archived locally in ~/FRO/Parsers/Effelsberg/Surveys/Holmberg1/Holmberg1 Raw/ (archive Holmberg1_Ef_data.tag.gz, ~3.4 GB uncompressed, plus the images Holmberg1.png/Holmberg1_zoom.png and the session's observation.log). The extracted MBFITS tree lives in Surveys/Holmberg1/Holmberg1_MBFITS_Profiles/."),
("H","6. Operating Procedure"),
("B","Step 1 — Extract the MBFITS archive."),
("P","Decompress Holmberg1_Ef_data.tag.gz into the Holmberg1_MBFITS_Profiles/ folder (the archive is a standard tar.gz despite the .tag.gz extension)."),
("B","Step 2 — Run the parser."),
("P","env -u LD_LIBRARY_PATH python3 ~/FRO/Parsers/Effelsberg/effelsberg_to_fro.py [MBFITS_DIR] [OUTPUT_DIR]. With no arguments it uses the Holmberg1 default paths. One .fro file is written per baseband into Holmberg1_FRO/."),
("B","Step 3 — Verify and archive."),
("P","Check the summary printed by read_fro_summary() (spectrum count, RA/DEC/GLON/GLAT range, non-tracking row count). Run the fro-backup alias to mirror ~/FRO/ into ~/FRO_Backup/."),
("H","7. Parser: effelsberg_to_fro.py"),
("P","Location: ~/FRO/Parsers/Effelsberg/. One .fro file per baseband (a single baseband in this scan), all subscans merged into the same file. Observation metadata (OBJECT, TELESCOP, SCANTYPE, SCANMODE, PROJID) is read from the primary header (HDU0) of GROUPING.fits inside the observation directory — not from a hypothetical external file, since in the archive the top-level summary .fits is actually a relative symlink to this same GROUPING.fits. facility='Effelsberg 100m RT' (from TELESCOP), target='HOI' (from OBJECT). polarization='pol1'/'pol2' for the two NUSEFEED channels. flag_cal = raw DATAPAR PHASE value (1.0 or 2.0). flag_track = 1 when the (LONGOFF, LATOFF) distance from the nominal position is below 1.0° (ON, constant ON_OFF_OFFSET_THRESHOLD_DEG in the script), 0 otherwise (OFF). Site coordinates and temperature-calibration notes are consistent with ebhis_to_fro.py in the same folder. Nothing is hardcoded: MBFITS path from argv or default, all frequency/polarization/pointing parameters read from the FITS headers."),
("H","8. ezRA Compatibility"),
("P","The internal HDF5 layout of format v17 (/Observatory, /Spectra, /Time, /Pointing) is the same one expected by fro_to_ezra.py, so .fro files produced by this parser convert to ezRA/ezCon.py's text format without changes to the conversion script."),
("H","9. References"),
("P","FRO format specification: fro_format_v17.py (~/FRO/FRO_System/)."),
("P","MBFITS format: same standard described in APEX_MBFITS_Data_Structure_EN_v1.1.docx (~/FRO/Parsers/APEX/), shared between APEX and Effelsberg."),
("P","MPIfR Effelsberg 100m: Max-Planck-Institut für Radioastronomie, Bonn, Germany."),
("P","Note: no specific Holmberg I bibliographic citations were included to avoid reporting unverified references in this session; to be added if Francesco specifies a preferred DOI/paper."),
]

def find_styles(doc):
    """Trova gli stili reali usati nel template SRT per titolo, heading, corpo, bold."""
    title_style = heading_style = body_style = None
    for p in doc.paragraphs:
        sn = p.style.name
        if sn.lower() in ("title","titolo","titolo principale") and title_style is None:
            title_style = p.style
        if sn.lower().startswith("heading") or sn.lower().startswith("titolo "):
            if heading_style is None: heading_style = p.style
    return title_style, heading_style

def build(template_path, content, out_path):
    from docx.enum.text import WD_ALIGN_PARAGRAPH
    doc = Document(template_path)
    orig = doc.paragraphs
    # allineamenti reali di titolo e sottotitolo dal template
    title_align = orig[0].alignment if len(orig) > 0 else WD_ALIGN_PARAGRAPH.CENTER
    sub_align = orig[1].alignment if len(orig) > 1 else WD_ALIGN_PARAGRAPH.CENTER
    # svuota il corpo
    for para in list(doc.paragraphs):
        para._element.getparent().remove(para._element)
    for kind, text in content:
        para = doc.add_paragraph()
        if kind == "T":
            try: para.style = doc.styles["Title"]
            except Exception: pass
            para.alignment = title_align
            para.add_run(text)
        elif kind == "S":
            para.alignment = sub_align
            r = para.add_run(text); r.italic = True
        elif kind == "H":
            try: para.style = doc.styles["Heading 1"]
            except Exception: pass
            para.add_run(text)
        elif kind == "B":
            r = para.add_run(text); r.bold = True
        else:
            para.add_run(text)
    doc.save(out_path)
    print("saved", out_path)

build(SRT_IT, IT, OUT_IT)
build(SRT_EN, EN, OUT_EN)
