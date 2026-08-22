import copy
from docx import Document

SRT_IT = "/home/franz/FRO/Parsers/SRT/SRT_Data_Structure_IT.docx"
SRT_EN = "/home/franz/FRO/Parsers/SRT/SRT_Data_Structure_EN.docx"
OUT_IT = "/home/franz/FRO/Parsers/JCMT/Surveys/COHRS_30p00_0p00_CUBE_3T2_R2/COHRS_Data_Structure_IT_v1.1.docx"
OUT_EN = "/home/franz/FRO/Parsers/JCMT/Surveys/COHRS_30p00_0p00_CUBE_3T2_R2/COHRS_Data_Structure_EN_v1.1.docx"

# ---- contenuto: (tipo, testo) ----
# T=title, S=sub, H=heading, P=paragrafo, B=paragrafo tutto grassetto
IT = [
("T","Struttura Dati COHRS e Modalità Operativa"),
("S","FRO / HALO (Hydrogen Atomic Line Observatory) — Francesco Di Giovanni, Bolzano, Italia"),
("S","Versione documento: 1.0 — 31 luglio 2026"),
("H","1. Panoramica della Survey"),
("P","Il parser cohrs_to_fro.py estrae spettri singoli dai cubi FITS posizione-posizione-velocità di COHRS e li converte nel formato FRO HDF5 (Hierarchical Data Format 5) v17. COHRS (12CO(3-2) High-Resolution Survey of the Galactic Plane) è una survey ad ampia copertura del piano galattico settentrionale nella transizione rotazionale 12CO J=3-2 a 345.796 GHz, realizzata con il JCMT (James Clerk Maxwell Telescope) tramite l’array eterodino HARP. A differenza dei parser single-dish che leggono un file di osservazione, questo parser legge uno spettro alla volta da un cubo 3D, a coordinate galattiche scelte interattivamente dall’utente. Un visualizzatore complementare, view_cohrs.py, produce mappe di cielo e diagrammi l-v direttamente dal cubo. Parser e viewer sono stati sviluppati e validati sul tile G30 (regione W43, braccio di Scudo)."),
("H","2. Telescopio e Strumento"),
("P","Tipo: telescopio submillimetrico single-dish da 15 m. Sito: Maunakea, Hawaii, USA — quota ~4092 m s.l.m. Coordinate: lat 19.8228°, lon -155.4770°. Operatore: East Asian Observatory (già Joint Astronomy Centre). Strumento: HARP (Heterodyne Array Receiver Programme), array a 16 pixel sul piano focale operante nella banda 325-375 GHz (banda B). Beam FWHM (Full Width at Half Maximum, larghezza a metà altezza) a 345 GHz: ~14 arcsec nativo; i cubi COHRS DR2 sono grigliati a 16.6 arcsec. Backend: ACSIS (Auto-Correlation Spectral Imaging System). Modalità osservativa: mappatura on-the-fly (OTF) a intreccio (basket-weave). La survey combina osservazioni di più semestri JCMT in un mosaico continuo."),
("H","3. Formato del Cubo"),
("P","I prodotti dati COHRS DR2 sono cubi FITS 3D con assi (longitudine galattica, latitudine galattica, velocità radiale). Le parole chiave WCS (World Coordinate System) rilevanti sono lette dinamicamente dall’header: CTYPE1=GLON-TAN, CTYPE2=GLAT-TAN, CTYPE3=VRAD con CUNIT3=km/s. L’asse spettrale è quindi già una velocità radiale nel sistema LSRK (Local Standard of Rest, Kinematic) — nessuna conversione di frequenza è necessaria per leggerlo. BUNIT=K (temperatura d’antenna T_A*). RESTFRQ=345.796 GHz. I cubi sono divisi in tile da 0.5° in longitudine; il tile di riferimento è 303 x 603 x 788 pixel (GLon x GLat x velocità), ricampionato a 0.635 km/s per canale nell’intervallo -200/+300 km/s. Un HDU VARIANCE accompagna i dati primari. Nota: CRPIX1 può cadere molto fuori dal cubo (il WCS è ancorato a una longitudine di riferimento comune alla survey); astropy.wcs lo gestisce correttamente."),
("H","4. Dataset di Riferimento"),
("B","Perché questo target:"),
("P","Il tile G30 (COHRS_30p00_0p00_CUBE_3T2_R2.fit) copre GLon ~29.75-30.25°, GLat ±0.5°, una delle regioni più ricche del piano galattico settentrionale. Si trova vicino alla tangente del braccio di Scudo-Centauro e include W43, un complesso molecolare \"mini-starburst\", rendendolo un campo di validazione ideale dove una mappa momento-0 e un diagramma l-v mostrano struttura reale invece che rumore."),
("B","Cosa ci si aspettava:"),
("P","Forte emissione 12CO(3-2) concentrata sul piano galattico (b~0), con nubi molecolari a velocità LSRK positive (grosso modo +40/+120 km/s) corrispondenti al gas lungo la linea di vista attraverso il braccio di Scudo. Una mappa momento-0 dovrebbe rivelare il tessuto delle nubi molecolari lungo il piano."),
("B","Cosa è stato ottenuto:"),
("P","Viewer e parser validati end-to-end. La mappa momento-0 mostra la banda di emissione CO concentrata sul piano, con complessi molecolari brillanti (W_CO fino a ~130 K km/s) allineati lungo b~0. Uno spettro singolo estratto a (GLon=29.95, GLat=-0.05) mostra un chiaro picco CO vicino a +100 km/s (doppia componente a ~+95 e ~+100 km/s), baseline piatta e rumore ~±1 K — esattamente il gas del braccio di Scudo atteso. La conversione velocità->frequenza->velocità (il parser scrive un asse in frequenza in Hz, il viewer riconverte in velocità) chiude il cerchio senza artefatti: il picco appare dove deve stare."),
("P","Caratteristiche tecniche: cubo 303 x 603 x 788 (GLon x GLat x velocità), float32. Larghezza canale 0.635 km/s, intervallo velocità -200/+300 km/s LSRK. BUNIT=K (T_A*). Frequenza a riposo 345.796 GHz (12CO J=3-2). Zero NaN nel tile di riferimento; picco di brillanza ~34 K. Passo griglia 0.00167° (pixel ~6 arcsec)."),
("H","5. Accesso ai Dati"),
("P","COHRS DR2 (Complete Data Release; Park et al. 2023, ApJS 264, 16) è pubblicamente disponibile su CANFAR (Canadian Astronomy Data Centre), DOI 10.11570/22.0078. La directory dati contiene: RELEASE2_CUBE_REBIN/ (i cubi 3D ricampionati a 0.635 km/s, un FITS per tile di 0.5° in longitudine, ~1.1 GB ciascuno), RELEASE2_INTEG/ (mappe momento-0 di intensità integrata già pronte), RELEASE2_LV/ (diagrammi longitudine-velocità già pronti), più mappe combinate dell’intera survey e un README. Scaricare i singoli tile come FITS diretto. Verificare la dimensione del file prima di scaricare: ogni cubo ricampionato è ~1.1 GB."),
("H","6. Modalità Operativa"),
("P","Procedura completa dall’identificazione del tile all’output FRO e ai plot:"),
("B","Passo 1 — Identificare il tile."),
("P","Scegliere una longitudine galattica di interesse lungo il piano settentrionale (COHRS copre GLon 9.5-62.3°). Le regioni molecolari più ricche si trovano verso GLon ~24-40° (braccio di Scudo, W43, W49). I nomi dei tile codificano le coordinate centrali: COHRS_XXpXX_YpYY_CUBE_3T2_R2.fit."),
("B","Passo 2 — Scaricare il cubo."),
("P","Dalla directory dati CANFAR (DOI 10.11570/22.0078), aprire RELEASE2_CUBE_REBIN/, verificare la dimensione del file (~1.1 GB) e scaricare il FITS del tile scelto. Salvarlo in ~/FRO/Parsers/JCMT/Surveys/."),
("B","Passo 3 — Avviare il parser."),
("P","Fare doppio clic sul launcher COHRS to FRO sul desktop, oppure eseguire: env -u LD_LIBRARY_PATH python3 ~/FRO/Parsers/JCMT/cohrs_to_fro.py <cubo.fit>. Il parser ricava automaticamente la directory di lavoro dal nome del file cubo: ~/FRO/Parsers/JCMT/Surveys/<tag>/FRO/."),
("B","Passo 4 — Estrarre gli spettri."),
("P","Quando richiesto, inserire GLon e GLat (gradi) di ogni punto di interesse. Il parser estrae lo spettro al pixel più vicino, registra nelle Note le coordinate richieste vs quelle effettive del centro pixel (onestà dei metadati) e scrive un file .fro per punto. Ripetere interattivamente per quanti punti si desidera."),
("B","Passo 5 — Generare i plot."),
("P","A fine estrazione il parser offre di eseguire view_cohrs.py sullo stesso cubo. Il menu del viewer offre: plot degli spettri estratti, diagramma l-v (a una GLat scelta o integrato su b), mappa momento-0 di intensità integrata W(CO), mappa della temperatura d’antenna di picco, mappa della velocità del picco, o tutti. Colormap e contorni sono scelti a runtime. Tutti i PNG sono salvati in <tag>/Plots/ e aperti automaticamente con xdg-open."),
("B","Passo 6 — Verificare e archiviare."),
("P","Verificare che spettri e mappe siano corretti. Eseguire l’alias rsync fro-backup per rispecchiare ~/FRO/ in ~/FRO_Backup/. Conservare il FITS originale del cubo come archivio raw permanente."),
("H","7. Parser: cohrs_to_fro.py"),
("P","Posizione: ~/FRO/Parsers/JCMT/. Launcher: ~/FRO/Parsers/JCMT/cohrs_to_fro_launcher.sh. Viewer: ~/FRO/Parsers/JCMT/view_cohrs.py. Estrazione interattiva di punti da un cubo 3D: l’utente fornisce GLon/GLat, il parser estrae lo spettro del pixel più vicino. La facility è fissa a JCMT 15m / HARP (COHRS è una survey a singolo telescopio), con coordinate di Maunakea. L’asse di velocità (CTYPE3=VRAD, CUNIT3=km/s, LSRK) è convertito in un asse di frequenza per il gruppo Spectra di FRO via f = f0*(1 - v/c), f0 = RESTFRQ dall’header. L’epoca di osservazione è impostata al sentinella esplicito UNKNOWN_EPOCH: COHRS è un mosaico OTF che combina più semestri JCMT, quindi l’epoca di un singolo pixel non è definita (stessa gestione introdotta per HI4PI in fro_format_v17). polarization=\"unknown\", flag_track=-1 (sconosciuto), pointing_source=\"measured\". Niente è hardcoded: percorso del cubo da argv o input(), directory di lavoro derivata dal nome del cubo, tutti i parametri WCS letti dall’header. Un file .fro per punto estratto."),
("H","8. Note Scientifiche"),
("P","La mappa momento-0 è presentata come intensità integrata W(CO) in K km/s, NON come densità di colonna H2: convertire in densità di colonna richiede un fattore di conversione CO-H2 assunto (X_CO), un’assunzione forte e incerta. W(CO) è la quantità onesta, misurata direttamente. Nel tile di riferimento G30 l’emissione CO è confinata a velocità LSRK positive (gas del braccio di Scudo); lo spettro a (29.95, -0.05) ha il picco vicino a +100 km/s con una doppia componente risolta. COHRS DR2 registra la temperatura d’antenna T_A*; un fattore di efficienza main-beam (~0.61, dal README della survey) converte in T_mb ma non è applicato dal parser. Per la mappa della velocità del picco, i pixel il cui picco è sotto una piccola soglia sono mascherati, altrimenti le linee di vista dominate dal rumore dipingono velocità casuali. Questa è la prima riga molecolare (CO 345 GHz) integrata nella pipeline FRO, a dimostrazione che FRO non è limitato alla riga HI ma è aperto a qualsiasi riga spettrale."),
("H","9. Campionamento Nyquist"),
("P","Grigliatura COHRS DR2: 16.6 arcsec (beam HARP nativo ~14 arcsec a 345 GHz). Il cubo è già campionato al Nyquist dalla survey. Per estrarre spettri indipendenti, punti distanziati di ~1 pixel (griglia ~6 arcsec) sono fortemente sovracampionati; beam indipendenti sono separati di ~14-17 arcsec. Nel mappare la struttura, scegliere punti di estrazione non più vicini della dimensione del beam per evitare spettri ridondanti e correlati."),
("H","10. Riferimenti"),
("P","Park, G., Currie, M. J., Thomas, H. S., et al. 2023, ApJS 264, 16 — COHRS Complete Data Release (DR2), il dataset processato da questo parser."),
("P","Dempsey, J. T., Thomas, H. S., Currie, M. J. 2013, ApJS 209, 8 — primo rilascio dati COHRS (DR1)."),
("P","Buckle, J. V., Hills, R. E., Smith, H., et al. 2009, MNRAS 399, 1026 — array eterodino HARP/ACSIS e backend al JCMT."),
("P","COHRS DR2 su CANFAR: DOI 10.11570/22.0078"),
("P","Specifica formato FRO: fro_format_v17.py (~/FRO/FRO_System/)"),
]

EN = [
("T","COHRS Data Structure and Operating Procedure"),
("S","FRO / HALO (Hydrogen Atomic Line Observatory) — Francesco Di Giovanni, Bolzano, Italy"),
("S","Document version: 1.0 — 31 July 2026"),
("H","1. Survey Overview"),
("P","The parser cohrs_to_fro.py extracts individual spectra from COHRS position-position-velocity FITS cubes and converts them into FRO HDF5 (Hierarchical Data Format 5) v17 format. COHRS (12CO(3-2) High-Resolution Survey of the Galactic Plane) is a large-area survey of the northern Galactic plane in the 12CO J=3-2 rotational transition at 345.796 GHz, carried out with the JCMT (James Clerk Maxwell Telescope) using the HARP heterodyne array. Unlike single-dish parsers that read one observation file, this parser reads one spectrum at a time from a 3D cube, at galactic coordinates chosen interactively by the user. A companion viewer, view_cohrs.py, produces sky maps and l-v diagrams directly from the cube. The parser and viewer were developed and validated on the G30 tile (W43 region, Scutum arm)."),
("H","2. Telescope and Instrument"),
("P","Type: 15-m submillimetre single-dish telescope. Site: Maunakea, Hawaii, USA — altitude ~4092 m a.s.l. Coordinates: lat 19.8228°, lon -155.4770°. Operator: East Asian Observatory (formerly Joint Astronomy Centre). Instrument: HARP (Heterodyne Array Receiver Programme), a 16-pixel focal-plane array operating in the 325-375 GHz band (B-band). Beam FWHM (Full Width at Half Maximum) at 345 GHz: ~14 arcsec native; COHRS DR2 cubes are gridded to 16.6 arcsec. Backend: ACSIS (Auto-Correlation Spectral Imaging System). Observing mode: basket-weave on-the-fly (OTF) mapping. The survey combines observations from multiple JCMT semesters into a seamless mosaic."),
("H","3. Cube Format"),
("P","COHRS DR2 data products are 3D FITS cubes with axes (Galactic longitude, Galactic latitude, radial velocity). The relevant WCS (World Coordinate System) keywords are read dynamically from the header: CTYPE1=GLON-TAN, CTYPE2=GLAT-TAN, CTYPE3=VRAD with CUNIT3=km/s. The spectral axis is therefore already a radial velocity in the LSRK (Local Standard of Rest, Kinematic) frame — no frequency conversion is needed to read it. BUNIT=K (antenna temperature T_A*). RESTFRQ=345.796 GHz. Cubes are divided into tiles of 0.5° in longitude; the reference tile is 303 x 603 x 788 pixels (GLon x GLat x velocity), rebinned to 0.635 km/s per channel over the range -200 to +300 km/s. A VARIANCE HDU accompanies the primary data. Note: CRPIX1 may fall far outside the cube (the WCS is anchored to a common survey reference longitude); astropy.wcs handles this correctly."),
("H","4. Reference Dataset"),
("B","Why this target:"),
("P","The G30 tile (COHRS_30p00_0p00_CUBE_3T2_R2.fit) covers GLon ~29.75-30.25°, GLat ±0.5°, one of the richest regions of the northern Galactic plane. It lies near the tangent of the Scutum-Centaurus arm and includes W43, a \"mini-starburst\" molecular complex, making it an ideal validation field where a moment-0 map and an l-v diagram show real structure rather than noise."),
("B","What we expected:"),
("P","Strong 12CO(3-2) emission concentrated on the Galactic plane (b~0), with molecular clouds at positive LSRK velocities (roughly +40 to +120 km/s) corresponding to gas along the line of sight through the Scutum arm. A moment-0 map should reveal the molecular cloud fabric along the plane."),
("B","What we obtained:"),
("P","Viewer and parser validated end-to-end. The moment-0 map shows the band of CO emission concentrated on the plane, with bright molecular complexes (W_CO up to ~130 K km/s) aligned along b~0. A single spectrum extracted at (GLon=29.95, GLat=-0.05) shows a clear CO peak near +100 km/s (a double component at ~+95 and ~+100 km/s), a flat baseline and ~±1 K noise — exactly the Scutum-arm gas expected. The velocity->frequency->velocity round trip (parser writes a frequency axis in Hz, viewer reconverts to velocity) closes without artefacts: the peak appears where it should."),
("P","Technical characteristics: cube 303 x 603 x 788 (GLon x GLat x velocity), float32. Channel width 0.635 km/s, velocity range -200 to +300 km/s LSRK. BUNIT=K (T_A*). Rest frequency 345.796 GHz (12CO J=3-2). Zero NaN in the reference tile; peak brightness ~34 K. Grid step 0.00167° (~6 arcsec pixels)."),
("H","5. Data Access"),
("P","COHRS DR2 (Complete Data Release; Park et al. 2023, ApJS 264, 16) is publicly available on CANFAR (Canadian Astronomy Data Centre), DOI 10.11570/22.0078. The data directory contains: RELEASE2_CUBE_REBIN/ (the 3D cubes rebinned to 0.635 km/s, one FITS per 0.5° longitude tile, ~1.1 GB each), RELEASE2_INTEG/ (ready-made moment-0 integrated-intensity maps), RELEASE2_LV/ (ready-made longitude-velocity diagrams), plus full-survey combined maps and a README. Download individual cube tiles directly as FITS. Verify file size before downloading: each rebinned cube is ~1.1 GB."),
("H","6. Operating Procedure"),
("P","Complete workflow from tile identification to FRO output and plots:"),
("B","Step 1 — Identify the tile."),
("P","Choose a Galactic longitude of interest along the northern plane (COHRS covers GLon 9.5-62.3°). Richer molecular regions lie toward GLon ~24-40° (Scutum arm, W43, W49). Tile names encode the central coordinates: COHRS_XXpXX_YpYY_CUBE_3T2_R2.fit."),
("B","Step 2 — Download the cube."),
("P","From the CANFAR data directory (DOI 10.11570/22.0078), open RELEASE2_CUBE_REBIN/, check the file size (~1.1 GB), and download the chosen tile FITS. Save it in ~/FRO/Parsers/JCMT/Surveys/."),
("B","Step 3 — Launch the parser."),
("P","Double-click the COHRS to FRO launcher on the desktop, or run: env -u LD_LIBRARY_PATH python3 ~/FRO/Parsers/JCMT/cohrs_to_fro.py <cube.fit>. The parser derives its working directory automatically from the cube file name: ~/FRO/Parsers/JCMT/Surveys/<tag>/FRO/."),
("B","Step 4 — Extract spectra."),
("P","When prompted, enter GLon and GLat (deg) of each point of interest. The parser extracts the spectrum at the nearest pixel, records requested vs actual pixel-centre coordinates in the Notes (metadata honesty), and writes one .fro file per point. Repeat interactively for as many points as wanted."),
("B","Step 5 — Generate plots."),
("P","At the end of extraction the parser offers to run view_cohrs.py on the same cube. The viewer menu offers: extracted-spectra plots, l-v diagram (at a chosen GLat or integrated over b), moment-0 integrated-intensity map W(CO), peak antenna-temperature map, velocity-of-peak map, or all of the above. Colormap and contour options are chosen at runtime. All PNGs are saved in <tag>/Plots/ and opened automatically with xdg-open."),
("B","Step 6 — Verify and archive."),
("P","Check that the spectra and maps look correct. Run the fro-backup rsync alias to mirror ~/FRO/ into ~/FRO_Backup/. Keep the original cube FITS as the permanent raw archive."),
("H","7. Parser: cohrs_to_fro.py"),
("P","Location: ~/FRO/Parsers/JCMT/. Launcher: ~/FRO/Parsers/JCMT/cohrs_to_fro_launcher.sh. Viewer: ~/FRO/Parsers/JCMT/view_cohrs.py. Interactive point extraction from a 3D cube: the user supplies GLon/GLat, the parser extracts the nearest-pixel spectrum. Facility is fixed to JCMT 15m / HARP (COHRS is a single-telescope survey), with Maunakea coordinates. The velocity axis (CTYPE3=VRAD, CUNIT3=km/s, LSRK) is converted to a frequency axis for the FRO Spectra group via f = f0*(1 - v/c), f0 = RESTFRQ from the header. Observation epoch is set to the explicit UNKNOWN_EPOCH sentinel: COHRS is an OTF mosaic combining several JCMT semesters, so the epoch of any single pixel is not defined (the same handling introduced for HI4PI in fro_format_v17). polarization=\"unknown\", flag_track=-1 (unknown), pointing_source=\"measured\". Nothing is hardcoded: cube path from argv or input(), working directory derived from the cube name, all WCS parameters read from the header. One .fro file per extracted point."),
("H","8. Scientific Notes"),
("P","The moment-0 map is presented as integrated intensity W(CO) in K km/s, NOT as an H2 column density: converting to column density requires an assumed CO-to-H2 conversion factor (X_CO), which is a strong and uncertain assumption. W(CO) is the honest, directly measured quantity. In the G30 reference tile, CO emission is confined to positive LSRK velocities (Scutum-arm gas); the spectrum at (29.95, -0.05) peaks near +100 km/s with a resolved double component. COHRS DR2 records antenna temperature T_A*; a main-beam efficiency factor (~0.61, from the survey README) converts to T_mb but is not applied by the parser. For the velocity-of-peak map, pixels whose peak lies below a small threshold are masked, otherwise noise-dominated lines of sight paint random velocities. This is the first molecular line (CO 345 GHz) integrated into the FRO pipeline, demonstrating that FRO is not limited to the HI line but is open to any spectral line."),
("H","9. Nyquist Sampling"),
("P","COHRS DR2 gridding: 16.6 arcsec (native HARP beam ~14 arcsec at 345 GHz). The cube is already Nyquist-sampled by the survey. For extracting independent spectra, points spaced by ~1 pixel (~6 arcsec grid) are heavily oversampled; independent beams are separated by ~14-17 arcsec. When mapping structure, choose extraction points no closer than roughly the beam size to avoid redundant, correlated spectra."),
("H","10. References"),
("P","Park, G., Currie, M. J., Thomas, H. S., et al. 2023, ApJS 264, 16 — COHRS Complete Data Release (DR2), the dataset processed by this parser."),
("P","Dempsey, J. T., Thomas, H. S., Currie, M. J. 2013, ApJS 209, 8 — COHRS first data release (DR1)."),
("P","Buckle, J. V., Hills, R. E., Smith, H., et al. 2009, MNRAS 399, 1026 — HARP/ACSIS heterodyne array and backend at the JCMT."),
("P","COHRS DR2 on CANFAR: DOI 10.11570/22.0078"),
("P","FRO format specification: fro_format_v17.py (~/FRO/FRO_System/)"),
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