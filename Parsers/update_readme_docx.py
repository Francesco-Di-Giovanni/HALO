# update_readme_docx.py
# One-off script: inserts the "Convenzioni di Organizzazione" section into
# FRO_Parsers_README.docx (IT), and builds FRO_Parsers_README_EN.docx as a
# full English translation of the whole document (including the new section).
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy

from docx import Document
import copy

IT_PATH = '/home/franz/FRO/Parsers/FRO_Parsers_README.docx'
EN_PATH = '/home/franz/FRO/Parsers/FRO_Parsers_README_EN.docx'

# ---------------------------------------------------------------------------
# New section content (IT), inserted before "Osservatori Supportati".
# kind: T=title (section header, like "Osservatori Supportati"),
#       B=bold subheading line, P=body paragraph, ''=blank spacer
# ---------------------------------------------------------------------------
NEW_SECTION_IT = [
    ('T', 'Convenzioni di Organizzazione del Progetto'),
    ('P', "I parser di questa directory sono organizzati per dataset/survey, non per telescopio. "
          "Questa sezione spiega il criterio e come si traduce in struttura di directory."),
    ('', ''),
    ('B', 'Telescopio, survey e osservazione puntata: tre concetti distinti'),
    ('P', "Telescopio: lo strumento fisico (Effelsberg 100m, ATCA, Parkes, JCMT, APEX, ...). Definisce sito, "
          "coordinate, hardware e il formato dei dati grezzi (MBFITS, SDFITS, RPFITS, FITS generico, ...). "
          "Una cartella Parsers/<Telescopio>/ raggruppa tutti gli script che leggono quel formato/hardware."),
    ('P', "Survey: un programma di osservazione sistematico e pianificato del cielo (o di un campione definito "
          "di oggetti), con una propria identità di dataset — spesso un data release pubblico, una propria "
          "pipeline di riduzione, un proprio nome — indipendentemente dal telescopio con cui è stata realizzata. "
          "Esempio: EBHIS è una survey HI di tutto il cielo boreale realizzata con Effelsberg, ma è un prodotto "
          "dati indipendente (profili ASCII già ridotti, scaricati dall'API del server di Bonn); ha quindi un "
          "parser proprio, ebhis_to_fro.py, distinto da quello che legge le scansioni MBFITS grezze dello stesso "
          "telescopio."),
    ('P', "Osservazione puntata: un singolo puntamento su un oggetto celeste specifico, non un programma "
          "sistematico. Esempio: Holmberg1 è una singola scansione MBFITS di Effelsberg su una galassia nana; "
          "NGC253 è una singola osservazione ATCA su una galassia. Un'osservazione puntata normalmente riusa il "
          "parser generico del formato/strumento (effelsberg_to_fro.py, atca_to_fro.py) invece di richiederne uno "
          "dedicato, perché non è un prodotto dati indipendente ma una singola acquisizione nello stesso formato "
          "già gestito."),
    ('', ''),
    ('B', 'Criterio per scrivere un parser dedicato'),
    ('P', "Si scrive un nuovo script parser quando cambia il formato o la provenienza dei dati grezzi (es. profili "
          "ASCII scaricati via web vs scansioni binarie MBFITS dallo stesso telescopio), non semplicemente perché "
          "cambia il target osservativo. Esempi: GASS (Parkes) ha parkes_gass_to_fro.py perché legge profili GASS "
          "già ridotti, non scansioni RPFITS grezze; COHRS (JCMT) ha cohrs_to_fro.py perché legge cubi FITS 3D "
          "pubblici della survey, non registrazioni ACSIS grezze."),
    ('', ''),
    ('B', 'Regole di nomenclatura delle directory'),
    ('P', "Parsers/<Telescopio>/ — cartella di primo livello per strumento fisico (APEX, ATCA, Effelsberg, FAST, "
          "GBT, JCMT, Parkes, SALSA, SRT, ...). Contiene gli script parser, i launcher .sh, l'icona .svg e la "
          "docx Data Structure (DS) generica dello strumento/formato."),
    ('P', "Parsers/<Telescopio>/Surveys/ — cartella dati, organizzata per survey o osservazione puntata:"),
    ('P', "— Survey sistematica: <NomeSurvey>_Survey_<tag>/ (es. EBHIS_Survey_G070-100_B-20+20, "
          "Parkes_Survey_GASS), con sottocartelle <TEL>_Profiles_<tag>/ (dati grezzi/scaricati), "
          "<TEL>_FRO_<tag>/ (output .fro) e <TEL>_Plots_<tag>/ (plot)."),
    ('P', "— Osservazione puntata: <NomeOggetto>/ (es. NGC253, Holmberg1) — struttura più piatta: cartella dati "
          "grezzi (spesso \"<NomeOggetto> Raw/\" per un archivio pesante da preservare), output .fro, e Plots/."),
    ('P', "Docx: la Data Structure (DS, panoramica formato e procedura) vive in Parsers/<Telescopio>/ se descrive "
          "un formato/parser generico, oppure accanto al dataset se è specifica di una survey/osservazione (es. "
          "COHRS_Data_Structure in JCMT/, Holmberg1_Data_Structure in Effelsberg/). Il Plot Description (PD) vive "
          "dentro la cartella Plots/ del dataset specifico (es. ATCA_Plots_Description_NGC253 in "
          "Surveys/NGC253/Plots/)."),
    ('P', "Dati raw: conservati integralmente nella cartella del dataset, mai sovrascritti — tipicamente in una "
          "sottocartella \"Raw\" per le singole osservazioni puntate con archivio pesante, oppure nella cartella "
          "*_Profiles_<tag>/ per le survey a griglia."),
    ('P', "Plot: PNG a 300 DPI per le mappe 2D, SVG per gli spettri 1D, HTML interattivi Plotly quando presenti — "
          "sempre in una sottocartella Plots/ del dataset."),
    ('', ''),
    ('B', 'Esempi dal progetto'),
    ('P', "EBHIS: survey Effelsberg con parser proprio (ebhis_to_fro.py), perché è un dataset indipendente "
          "(profili ASCII pubblici, non scansioni MBFITS)."),
    ('P', "Holmberg1: osservazione puntata Effelsberg, riusa la logica generica MBFITS (effelsberg_to_fro.py) "
          "perché non è una survey a sé stante."),
    ('P', "GASS: survey Parkes con parser proprio (parkes_gass_to_fro.py), stesso criterio di EBHIS."),
    ('P', "NGC253: osservazione puntata ATCA — cartella Surveys/NGC253/ piatta, senza la suddivisione "
          "Profiles/FRO/Plots annidata tipica delle survey a griglia."),
    ('P', "COHRS: survey JCMT con parser proprio (cohrs_to_fro.py), che legge cubi 3D pubblici già ridotti — "
          "distinto da un ipotetico parser di dati grezzi ACSIS dello stesso telescopio."),
    ('', ''),
]

NEW_SECTION_EN = [
    ('T', 'Project Organization Conventions'),
    ('P', "The parsers in this directory are organized by dataset/survey, not by telescope. This section "
          "explains the criterion and how it translates into directory structure."),
    ('', ''),
    ('B', 'Telescope, survey and pointed observation: three distinct concepts'),
    ('P', "Telescope: the physical instrument (Effelsberg 100m, ATCA, Parkes, JCMT, APEX, ...). It defines the "
          "site, coordinates, hardware and raw data format (MBFITS, SDFITS, RPFITS, generic FITS, ...). One "
          "Parsers/<Telescope>/ folder groups all scripts that read that format/hardware."),
    ('P', "Survey: a systematic, planned sky observation program (or of a defined object sample), with its own "
          "dataset identity — often a public data release, its own reduction pipeline, its own name — "
          "independent of the telescope it was carried out with. Example: EBHIS is a northern all-sky HI survey "
          "carried out with Effelsberg, but it is an independent data product (already-reduced ASCII profiles "
          "downloaded from the Bonn server API); it therefore has its own parser, ebhis_to_fro.py, distinct from "
          "the one reading raw MBFITS scans from the same telescope."),
    ('P', "Pointed observation: a single pointing at one specific celestial object, not a systematic program. "
          "Example: Holmberg1 is a single Effelsberg MBFITS scan of a dwarf galaxy; NGC253 is a single ATCA "
          "observation of a galaxy. A pointed observation normally reuses the generic format/instrument parser "
          "(effelsberg_to_fro.py, atca_to_fro.py) instead of requiring a dedicated one, because it is not an "
          "independent data product but a single acquisition in an already-handled format."),
    ('', ''),
    ('B', 'Criterion for writing a dedicated parser'),
    ('P', "A new parser script is written when the format or provenance of the raw data changes (e.g. ASCII "
          "profiles downloaded via the web vs binary MBFITS scans from the same telescope), not simply because "
          "the observing target changes. Examples: GASS (Parkes) has parkes_gass_to_fro.py because it reads "
          "already-reduced GASS profiles, not raw RPFITS scans; COHRS (JCMT) has cohrs_to_fro.py because it "
          "reads public 3D FITS cubes from the survey, not raw ACSIS recordings."),
    ('', ''),
    ('B', 'Directory naming rules'),
    ('P', "Parsers/<Telescope>/ — top-level folder per physical instrument (APEX, ATCA, Effelsberg, FAST, GBT, "
          "JCMT, Parkes, SALSA, SRT, ...). Contains the parser scripts, the .sh launchers, the .svg icon, and "
          "the generic instrument/format Data Structure (DS) docx."),
    ('P', "Parsers/<Telescope>/Surveys/ — data folder, organized by survey or pointed observation:"),
    ('P', "— Systematic survey: <SurveyName>_Survey_<tag>/ (e.g. EBHIS_Survey_G070-100_B-20+20, "
          "Parkes_Survey_GASS), with subfolders <TEL>_Profiles_<tag>/ (raw/downloaded data), "
          "<TEL>_FRO_<tag>/ (.fro output) and <TEL>_Plots_<tag>/ (plots)."),
    ('P', "— Pointed observation: <ObjectName>/ (e.g. NGC253, Holmberg1) — flatter structure: a raw-data folder "
          "(often \"<ObjectName> Raw/\" for a heavy archive to preserve), .fro output, and Plots/."),
    ('P', "Docx: the Data Structure (DS, format and procedure overview) lives in Parsers/<Telescope>/ if it "
          "describes a generic format/parser, or alongside the dataset if it is survey/observation-specific "
          "(e.g. COHRS_Data_Structure in JCMT/, Holmberg1_Data_Structure in Effelsberg/). The Plot Description "
          "(PD) lives inside the specific dataset's Plots/ folder (e.g. ATCA_Plots_Description_NGC253 in "
          "Surveys/NGC253/Plots/)."),
    ('P', "Raw data: kept in full inside the dataset folder, never silently overwritten — typically in a \"Raw\" "
          "subfolder for individual pointed observations with a heavy archive, or in the *_Profiles_<tag>/ "
          "folder for grid surveys."),
    ('P', "Plots: 300 DPI PNG for 2D maps, SVG for 1D spectra, interactive Plotly HTML where present — always "
          "in a Plots/ subfolder of the dataset."),
    ('', ''),
    ('B', 'Examples from the project'),
    ('P', "EBHIS: an Effelsberg survey with its own parser (ebhis_to_fro.py), because it is an independent "
          "dataset (public ASCII profiles, not MBFITS scans)."),
    ('P', "Holmberg1: an Effelsberg pointed observation, reusing the generic MBFITS logic "
          "(effelsberg_to_fro.py) because it is not a survey in its own right."),
    ('P', "GASS: a Parkes survey with its own parser (parkes_gass_to_fro.py), same criterion as EBHIS."),
    ('P', "NGC253: an ATCA pointed observation — a flat Surveys/NGC253/ folder, without the nested "
          "Profiles/FRO/Plots split typical of grid surveys."),
    ('P', "COHRS: a JCMT survey with its own parser (cohrs_to_fro.py), reading public, already-reduced 3D "
          "cubes — distinct from a hypothetical raw-ACSIS-data parser for the same telescope."),
    ('', ''),
]

# ---------------------------------------------------------------------------
# Full English translation of the existing per-observatory content, in the
# same order as the IT document (title, bold subheading, body paragraph).
# ---------------------------------------------------------------------------
OBSERVATORIES_EN = [
    ('APEX / MBFITS  —  Atacama Pathfinder EXperiment',
     "ESO/MPIfR/OSO 12m radio telescope on the Chilean Altiplano at 5100m a.s.l. The parser mbfits_to_fro.py "
     "converts MBFITS (Multi-Beam FITS) format data, the standard European format for single-dish radio "
     "telescopes. Reference dataset: 12CO (345 GHz) observation from project APEX-48698."),
    ('ATCA  —  Australia Telescope Compact Array',
     "Australian 6-antenna array operated by ATNF (CSIRO), located at Narrabri, New South Wales. Covers the "
     "HI band at 1420 MHz with high angular resolution. Raw data is in RPFITS format (ATNF-specific binary, "
     "interferometric visibility data). Reference dataset: NGC 253, project C1025, ATOA archive. Parser: "
     "placeholder — calibration requires CASA."),
    ('EBHIS  —  Effelsberg-Bonn HI Survey',
     "All-sky northern-hemisphere HI survey carried out with the Effelsberg 100m radio telescope (MPIfR, "
     "Germany). Together with the Australian GASS it forms the HI4PI catalogue. The parser ebhis_to_fro.py "
     "downloads and converts ASCII profiles from the Bonn server API (astro.uni-bonn.de). Reference dataset: "
     "GLon 70-100 deg, GLat +/-20 deg, 2.5 deg step."),
    ('FAST — raw  —  Five-hundred-meter Aperture Spherical Telescope',
     "The world's largest single-dish radio telescope, 500m, located in Guizhou, China (NAOC). The parser "
     "fast_fits_to_fro.py converts raw OTF (On-The-Fly) data in FITS format, integrating pointing coordinates "
     "from the KY table (.xlsx) via FAST_KY_to_radec.py (Dr. Jing Wang, NAOC). Reference dataset: M33, beam "
     "M01, 2021-07-30."),
    ('FAST — FEASTS  —  Furthest gas in FAST, Extending Armlength To Study HI',
     "HI survey of 55 Local Volume galaxies with FAST, reduced and calibrated data (Jy/beam) from the HiFAST "
     "v1.2 pipeline. The parser feasts_cube_to_fro.py extracts spectra from 3D FITS cubes (RA x Dec x VRAD) "
     "at interactively chosen positions. Reference dataset: NGC 628, 2022-11-04. Public data: "
     "disk.pku.edu.cn/link/AA5A6A26820B154A2F835BE8BE0FD906FD (Wang et al. 2025)."),
    ('GBT  —  Green Bank Telescope',
     "The world's largest fully steerable radio telescope, 100m, in West Virginia, USA (NRAO). The parser "
     "gbt_sdfits_to_fro.py converts SDFITS (Single Dish FITS) format data. Reference dataset: UGC 11891 (HI, "
     "detection at ~1418 MHz) and Mrk 1 (22 GHz H2O maser), NASA FITS Registry public archive."),
    ('HI4PI  —  HI 4pi survey',
     "All-sky HI survey combining EBHIS (Effelsberg, northern hemisphere) and GASS (Parkes, southern "
     "hemisphere). Angular resolution 16.2 arcmin, distributed as 3D FITS cubes (CAR projection). The parser "
     "hi4pi_to_fro.py extracts spectra from position-position-velocity cubes. Reference dataset: tiles CAR_D05 "
     "and CAR_E18. Public data: cdsarc.cds.unistra.fr."),
    ('LAB  —  Leiden/Argentine/Bonn survey',
     "Historic (2005) all-sky HI survey, 0.6 deg beam, ~1 km/s spectral resolution. The parser lab_to_fro.py "
     "downloads and converts ASCII profiles from the Bonn server API. Useful as a cross-check against "
     "EBHIS/HI4PI for FRO spectrum validation. Reference dataset: GLon 70-100 deg, GLat +/-20 deg, 2.5 deg "
     "step."),
    ('MeerKAT  —  MeerKAT Radio Telescope',
     "South African 64-antenna array, SKA precursor, located in the Karoo (South Africa). Excellent HI "
     "sensitivity, rich public archive (SARAO). Parser: placeholder — no sample file available at this time."),
    ('Medicina  —  Medicina Radio Astronomical Station',
     "INAF 32m radio telescope located near Bologna. Part of the Italian VLBI network together with SRT and "
     "Noto. Data accessible via the INAF archive and the EVN Data Archive (JIVE, Netherlands). Parser: "
     "placeholder — no sample file available at this time."),
    ('Noto  —  Noto Radio Astronomical Station',
     "INAF 32m radio telescope located in Sicily. Part of the Italian VLBI network together with SRT and "
     "Medicina. Data accessible via the INAF archive and the EVN Data Archive (JIVE, Netherlands). Parser: "
     "placeholder — no sample file available at this time."),
    ('Parkes  —  Parkes Radio Telescope — "The Dish"',
     "Australian 64m radio telescope operated by ATNF (CSIRO). Carried out the GASS survey (Galactic All-Sky "
     "Survey, part of HI4PI) and the HIPASS survey (HI Parkes All Sky Survey). The parser "
     "parkes_gass_to_fro.py converts GASS profiles. Reference dataset: acquisition in progress."),
    ('SALSA  —  Such A Lovely Small Antenna',
     "Network of small educational radio telescopes at the Onsala Space Observatory (Sweden), remotely "
     "operable via the web. Our primary instrument for the Galactic HI Survey at 1420 MHz. The parser "
     "salsa_fits_to_fro.py converts SALSA FITS observations. Reference dataset: Galactic survey GLon "
     "120-200 deg, GLat +/-20 deg, 5 deg step."),
    ('SRT  —  Sardinia Radio Telescope',
     "Italian 64m radio telescope operated by INAF at Pranu Sanguni, Sardinia. The parser srt_fits_to_fro.py "
     "converts data in fitszilla format (multi-extension FITS produced by SDT). Part of the Italian INAF VLBI "
     "network together with Medicina and Noto. Reference dataset: 3C286, L band, 2019 (restricted)."),
    ('THINGS  —  The HI Nearby Galaxy Survey',
     "High angular resolution HI survey of 34 nearby galaxies with the VLA (NRAO). 3D FITS cubes (RA x Dec x "
     "velocity). MPIA server currently offline — data not directly downloadable. Parser: placeholder."),
    ('VLA / VLASS  —  Very Large Array / VLA Sky Survey',
     "27-antenna array in New Mexico, USA (NRAO). VLASS is the 3 GHz continuum survey (S band, 2-4 GHz), "
     "~2.5 arcsec resolution, public Quick Look images via CADC (astroquery). The script view_vlass.py "
     "displays 2D FITS cutouts. Reference dataset: IGR J18249-3243 (FR II radio galaxy), 4 epochs 2017-2024."),
    ('WSRT  —  Westerbork Synthesis Radio Telescope',
     "Historic Dutch 14-antenna array (ASTRON), now partially replaced by APERTIF. Vast public HI archives. "
     "Parser: placeholder — no sample file available at this time."),
]

FOOTER_EN = 'Francesco Radio Observatory — FRO / HALO — 9 August 2026'


def insert_before(doc, ref_paragraph, kind, text):
    """Inserts a new paragraph of the given kind immediately before ref_paragraph."""
    new_p = ref_paragraph.insert_paragraph_before()
    if kind == 'T':
        new_p.style = doc.styles['Title']
        new_p.add_run(text)
    elif kind == 'B':
        new_p.style = doc.styles['Normal']
        r = new_p.add_run(text)
        r.bold = True
    elif kind == 'P':
        new_p.style = doc.styles['Normal']
        new_p.add_run(text)
    else:
        new_p.style = doc.styles['Normal']
    return new_p


# ---------------------------------------------------------------------------
# 1) Insert the new IT section into the existing README, before
#    "Osservatori Supportati".
# ---------------------------------------------------------------------------
doc_it = Document(IT_PATH)
target = None
for p in doc_it.paragraphs:
    if p.text.strip() == 'Osservatori Supportati':
        target = p
        break
if target is None:
    raise RuntimeError('Could not find "Osservatori Supportati" title paragraph')

for kind, text in NEW_SECTION_IT:
    insert_before(doc_it, target, kind, text)

doc_it.save(IT_PATH)
print('Updated', IT_PATH)

# ---------------------------------------------------------------------------
# 2) Build the full EN document, using the (now updated) IT docx as a style
#    template, with the main title, the new section, and every observatory
#    entry translated to English.
# ---------------------------------------------------------------------------
doc_en = Document(IT_PATH)
for para in list(doc_en.paragraphs):
    para._element.getparent().remove(para._element)

def add(kind, text):
    p = doc_en.add_paragraph()
    if kind == 'T':
        p.style = doc_en.styles['Title']
        p.add_run(text)
    elif kind == 'B':
        p.style = doc_en.styles['Normal']
        r = p.add_run(text)
        r.bold = True
    elif kind == 'I':
        p.style = doc_en.styles['Normal']
        r = p.add_run(text)
        r.italic = True
    elif kind == 'P':
        p.style = doc_en.styles['Normal']
        p.add_run(text)
    else:
        p.style = doc_en.styles['Normal']

add('T', 'FRO Parsers')
for kind, text in NEW_SECTION_EN:
    add(kind, text)
add('T', 'Supported Observatories')
add('P', "This directory contains the parsers for importing and converting data from various radio "
         "telescopes into the FRO format. Below is a description of each observatory, in alphabetical order.")
add('', '')
for header, body in OBSERVATORIES_EN:
    add('B', header)
    add('P', body)
    add('', '')
add('I', FOOTER_EN)

doc_en.save(EN_PATH)
print('Created', EN_PATH)
