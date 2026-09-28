#!/usr/bin/env python3
"""
classify_reads.py — VH/VL content deconvolution for one FASTQ shard.

Each read is aligned once against a combined index of designed VH-side and VL-side
half-references. A half is assigned only when its best hit beats the runner-up by
`--margin` matching bases over an alignment of at least `--min-mlen`; the read is then
classed by where the two halves come from.

Why content deconvolution rather than the aligner's own soft-clip / supplementary calls:
in a conserved-framework antibody library a chimera aligns end-to-end to one designed
gene with no soft clip, so a rescue-style caller scores it clean. The undercount is
library-specific and has been measured between ~9x and ~80x on other builds, so this
step must never be replaced by a soft-clip heuristic.

Emits ONE compact JSON per shard — class counts, per-gene and per-well counts, and read
length histograms — never per-read rows. A 500k-read shard would otherwise dominate the
work directory and the merge step.
"""
import argparse, gzip, json, os, sys
from collections import Counter, defaultdict

MARGIN, MIN_MLEN = 15, 120


def open_fastq(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def load_aligner(fa, best_n=5):
    import mappy
    return mappy.Aligner(fa, preset="map-ont", best_n=best_n)


def classify_seq(aln, seq, margin=MARGIN, min_mlen=MIN_MLEN):
    """Return (vh_gene|None, vl_gene|None, diagnostics)."""
    best_h, best_l = {}, {}
    for h in aln.map(seq):
        d = best_h if h.ctg[0] == "H" else best_l
        g = h.ctg[2:]
        if h.mlen > d.get(g, (0, None))[0]:
            d[g] = (h.mlen, h.strand)

    def pick(d):
        if not d:
            return None, 0, 0, None
        ranked = sorted(d.items(), key=lambda kv: -kv[1][0])
        g1, (m1, s1) = ranked[0]
        m2 = ranked[1][1][0] if len(ranked) > 1 else 0
        if m1 < min_mlen:
            return None, m1, m2, s1
        return (g1 if (m1 - m2) >= margin else None), m1, m2, s1

    vh, hm1, hm2, hs = pick(best_h)
    vl, lm1, lm2, ls = pick(best_l)
    return vh, vl, {"h_top": hm1, "h_run": hm2, "l_top": lm1, "l_run": lm2,
                    "strand_concordant": hs is not None and ls is not None and hs == ls}


def label(vh, vl, meta):
    if vh is None and vl is None:
        return "unmapped"
    if vh is None or vl is None:
        return "one_half"
    if vh == vl:
        return "clean_designed"
    a, b = meta.get(vh), meta.get(vl)
    if a is None or b is None:
        return "ambiguous"
    if a["order"] != b["order"]:
        return "cross_sublib"
    if a["well"] == b["well"]:
        return "within_well"
    return "cross_well"


PAIRED = ("clean_designed", "within_well", "cross_well", "cross_sublib")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--halves", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--sublibrary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--margin", type=int, default=MARGIN)
    ap.add_argument("--min-mlen", type=int, default=MIN_MLEN)
    ap.add_argument("--max-reads", type=int, default=0, help="0 = all reads (testing only)")
    a = ap.parse_args()

    meta = json.load(open(a.meta))
    aln = load_aligner(a.halves)
    cls = Counter()
    per_gene, per_well = Counter(), defaultdict(Counter)
    rlen_hist = Counter()                       # 25 bp bins, all reads
    rlen_by_cls = defaultdict(Counter)          # 25 bp bins, per class
    lenclass = defaultdict(Counter)             # single-molecule / fused counters
    n = 0
    with open_fastq(a.fastq) as fh:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            s = line.strip()
            if not s:
                continue
            n += 1
            if a.max_reads and n > a.max_reads:
                n -= 1
                break
            vh, vl, _ = classify_seq(aln, s, a.margin, a.min_mlen)
            c = label(vh, vl, meta)
            cls[c] += 1
            b = min(len(s), 2000) // 25 * 25
            rlen_hist[b] += 1
            rlen_by_cls[c][b] += 1
            lenclass[c]["n"] += 1
            if 600 <= len(s) <= 900:
                lenclass[c]["single_600_900"] += 1
            if len(s) >= 1200:
                lenclass[c]["fused_ge1200"] += 1
            lenclass[c]["len_sum"] += len(s)
            if c == "clean_designed":
                per_gene[vh] += 1
                per_well[meta[vh]["well"]][a.sublibrary] += 1
            if c in PAIRED and vh is not None:
                m = meta[vh]
                k = f"{m['order']}\t{m['well']}"
                per_well[k][c] += 1
                if c == "clean_designed" and len(s) >= 0.9 * m["full_len"]:
                    cls["_full_length_clean"] += 1
                if len(s) >= 0.9 * m["full_len"]:
                    cls["_full_length_assigned"] += 1

    out = {"sublibrary": a.sublibrary, "fastq": os.path.basename(a.fastq), "n_reads": n,
           "class_counts": dict(cls),
           "per_gene": dict(per_gene),
           "per_well": {k: dict(v) for k, v in per_well.items() if "\t" in k},
           "rlen_hist_25bp": dict(rlen_hist),
           "rlen_hist_25bp_by_cls": {k: dict(v) for k, v in rlen_by_cls.items()},
           "length_classes": {k: dict(v) for k, v in lenclass.items()},
           "params": {"margin": a.margin, "min_mlen": a.min_mlen}}
    json.dump(out, open(a.out, "w"))
    print(f"{a.sublibrary} {os.path.basename(a.fastq)}: {n} reads, "
          + ", ".join(f"{k} {v}" for k, v in sorted(cls.items()) if not k.startswith("_")))


if __name__ == "__main__":
    main()
