# fix_status_docx.py
# One-off script: fixes the "QR docx" header cell shading in
# FRO_Parsers_STATUS_v1.3.docx (was unfilled, must match the other header
# cells) and adds the Effelsberg/Holmberg1 row (alphabetically after EBHIS).
# FRO/HALO Project — Francesco Di Giovanni, Bolzano, Italy

from docx import Document
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
import copy

PATH = '/home/franz/FRO/Parsers/FRO_Parsers_STATUS_v1.3.docx'
HEADER_FILL = 'D3D1C7'

doc = Document(PATH)
table = doc.tables[0]

# --- 1) Fix "QR docx" header cell shading ---
header_row = table.rows[0]
qr_cell = None
for cell in header_row.cells:
    if cell.text.strip() == 'QR docx':
        qr_cell = cell
        break
if qr_cell is None:
    raise RuntimeError('"QR docx" header cell not found')

tcPr = qr_cell._tc.get_or_add_tcPr()
shd = tcPr.find(qn('w:shd'))
if shd is None:
    shd = OxmlElement('w:shd')
    tcPr.append(shd)
shd.set(qn('w:val'), 'clear')
shd.set(qn('w:color'), 'auto')
shd.set(qn('w:fill'), HEADER_FILL)
print('Fixed QR docx header cell fill ->', HEADER_FILL)

# --- 2) Add the Effelsberg/Holmberg1 row, after EBHIS ---
NEW_ROW = [
    'Effelsberg / Holmberg1',
    '✓ effelsberg_to_fro.py',
    '✗',
    '✓ Holmberg1 scan 6723',
    '✓',
    '✓ IT+EN',
    '✗',
    '✗',
    '✗',
]

ebhis_row = None
for row in table.rows:
    if row.cells[0].text.strip() == 'EBHIS':
        ebhis_row = row
        break
if ebhis_row is None:
    raise RuntimeError('EBHIS row not found')

# Clone the EBHIS row's XML (preserves borders/cell formatting), then
# overwrite its text and move it right after the EBHIS row.
new_tr = copy.deepcopy(ebhis_row._tr)
for cell_xml, text in zip(new_tr.findall(qn('w:tc')), NEW_ROW):
    # remove existing paragraphs, add a single fresh one with the new text
    for p in cell_xml.findall(qn('w:p')):
        cell_xml.remove(p)
    new_p = OxmlElement('w:p')
    new_r = OxmlElement('w:r')
    new_t = OxmlElement('w:t')
    new_t.text = text
    new_r.append(new_t)
    new_p.append(new_r)
    cell_xml.append(new_p)

ebhis_row._tr.addnext(new_tr)
print('Inserted Effelsberg / Holmberg1 row after EBHIS')

doc.save(PATH)
print('Saved', PATH)
