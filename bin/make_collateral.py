#!/usr/bin/env python3
"""
make_collateral.py — the one-page customer QC summary, in the MGF collateral format.

Deliberately numeric: five metrics and the read-length distribution, with a definitions
footer naming every denominator. Interpretation belongs in the full report, which a
scientist writes from these same outputs — this page is what goes to the customer with the
shipment, so it must be reproducible from the run alone.

Fields the run cannot know (account, order number, manufacture date) are highlighted as
placeholders rather than guessed.
"""
import json, os, sys
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_COLOR_INDEX
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn
from docx.oxml import OxmlElement

import argparse
AP = argparse.ArgumentParser(description="One-page MGF-style customer QC summary.")
AP.add_argument("--metrics-dir", default=".")
AP.add_argument("--plots-dir", default="plots")
AP.add_argument("--out", required=True)
AP.add_argument("--product", default="Direct-Print scFv Library (Ultra Long MGF)")
AP.add_argument("--account", default="")
AP.add_argument("--order-number", default="")
AP.add_argument("--order-items", default="")
AP.add_argument("--manufacture-date", default="")
AP.add_argument("--logo", default="")
A = AP.parse_args()
D = A.metrics_dir.rstrip("/") + "/"; F = A.plots_dir.rstrip("/") + "/"
OUT = A.out; MFG_DATE = A.manufacture_date or None
J = lambda n: json.load(open(D + n))
core, der, SC = J("metrics_core.json"), J("derived.json"), J("blockdel_scope.json")
ov, U = core["overall"], der["uniformity"]["ALL"]
chim = ov["within_well_pct_of_assigned"] + ov["cross_well_pct_of_assigned"] + ov["cross_sublib_pct_of_assigned"]
FL = SC.get("fl_matrix", {}).get("readlen_and_nogap", {}).get("both_primer")

doc = Document(); s = doc.sections[0]
s.page_width, s.page_height = Inches(8.5), Inches(11)
s.top_margin, s.bottom_margin, s.left_margin, s.right_margin = Inches(.6), Inches(.8), Inches(.9), Inches(.9)
st = doc.styles["Normal"]; st.font.name = "Arial"; st.font.size = Pt(10.5); st.paragraph_format.space_after = Pt(0)

def cellw(t, widths):
    t.autofit = False
    for gc, w in zip(t._tbl.tblGrid.findall(qn("w:gridCol")), widths): gc.set(qn("w:w"), str(int(w * 1440)))
    for r in t.rows:
        for c, w in zip(r.cells, widths): c.width = Inches(w)

def no_borders(t):
    tblPr = t._tbl.tblPr; b = OxmlElement("w:tblBorders")
    for e in ("top", "left", "bottom", "right", "insideH", "insideV"):
        x = OxmlElement(f"w:{e}"); x.set(qn("w:val"), "nil"); b.append(x)
    tblPr.append(b)

def rule(color="2BB39B", sz=8):
    p = doc.add_paragraph(); pPr = p._p.get_or_add_pPr(); bd = OxmlElement("w:pBdr"); bt = OxmlElement("w:bottom")
    bt.set(qn("w:val"), "single"); bt.set(qn("w:sz"), str(sz)); bt.set(qn("w:color"), color); bd.append(bt); pPr.append(bd)
    p.paragraph_format.space_after = Pt(10)

# letterhead (logo + address band, matching the existing customer QC summaries)
if A.logo and os.path.exists(A.logo):
    doc.add_paragraph().add_run().add_picture(A.logo, width=Inches(6.7))
    doc.add_paragraph().paragraph_format.space_after = Pt(6)
else:
    rule()

# order block
t = doc.add_table(rows=5, cols=2); no_borders(t); t.alignment = WD_TABLE_ALIGNMENT.CENTER; cellw(t, (2.6, 3.6))
ORD = sorted(core["per_sublib"])
rows = [("Product:", A.product), ("Account:", A.account or "[account]"), ("Order Number:", A.order_number or "[order number]"),
        ("Order Item Name:", A.order_items or ", ".join(ORD)), ("Manufacture Date:", MFG_DATE or "[confirm date]")]
for i, (k, v) in enumerate(rows):
    a = t.cell(i, 0).paragraphs[0]; a.alignment = WD_ALIGN_PARAGRAPH.RIGHT; r = a.add_run(k); r.bold = True
    b = t.cell(i, 1).paragraphs[0]; r = b.add_run("  " + v)
    if v.startswith("["): r.font.highlight_color = WD_COLOR_INDEX.YELLOW
    for c in (a, b): c.paragraph_format.space_after = Pt(5)
rule("BFBFBF", 6)

# metrics table
t = doc.add_table(rows=5, cols=2); t.style = "Table Grid"; cellw(t, (1.6, 5.1))
vals = [("Full Length Fragments", f"{FL:.1f}%" if FL is not None else "n/a"),
        ("95th/5th Percentile", f"{U['p95p5']:.2f}" if U.get("p95p5") else "n/a"),
        ("Dropouts", f"{U['dropout']} of {U['n_designed']:,} variants ({U['dropout_pct']:.2f}%)"),
        ("Pairing Chimera", f"{chim:.4f}%"), ("Read Length", None)]
for i, (k, v) in enumerate(vals):
    t.cell(i, 0).paragraphs[0].add_run(k)
    if v: t.cell(i, 1).paragraphs[0].add_run(v)
    else: t.cell(i, 1).paragraphs[0].add_run().add_picture(F + "collateral_readlen.png", width=Inches(4.9))
    for c in (t.cell(i, 0), t.cell(i, 1)):
        c.paragraphs[0].paragraph_format.space_before = Pt(3); c.paragraphs[0].paragraph_format.space_after = Pt(3)
rule("BFBFBF", 6)
p = doc.add_paragraph(); r = p.add_run(
    "Definitions. Full Length Fragments: reads ≥90% of designed length with no internal deletion ≥50 nt, among reads carrying both "
    "primer sites. Pairing Chimera: reads whose VH and VL halves come from different designed variants (within-well, cross-well and "
    "cross-sub-library combined), as a share of reads with both halves assigned. 95th/5th Percentile: ratio of per-variant read counts. "
    f"Full-depth Oxford Nanopore sequencing, {core['total_reads']/1e6:.1f} M reads across the three sub-libraries.")
r.font.size = Pt(8); r.font.color.rgb = RGBColor(0x59, 0x59, 0x59)
os.makedirs(os.path.dirname(os.path.abspath(OUT)) or ".", exist_ok=True)
doc.save(OUT); print("wrote", OUT)
