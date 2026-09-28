#!/usr/bin/env python3
"""
build_reference.py — turn the design reference CSV into the files the QC steps align against.

The reference CSV is the single design input to the pipeline. Everything downstream is
derived here, so a reference that passes this step is guaranteed self-consistent: the
half-references cannot drift out of step with the constructs, because both are cut from
the same `construct_dna` column using that row's own `vh_len_nt`.

Required columns: gene, order, well, vh_len_nt, construct_dna
Optional columns: framework, pool, well_rep2, length_nt, intermediate_len_nt, gc_pct

Outputs
  halves.fasta.gz      H_<gene> / L_<gene>, the classifier index
  constructs.fasta.gz  <gene>, the per-construct splice-alignment target
  gene_meta.json       {gene: {order, well, framework, vh_len, full_len, h_len, l_len}}
  reference_qc.json    what was checked, and anything that failed
"""
import argparse, csv, gzip, json, sys
from collections import Counter, defaultdict

DEFAULT_INT_LEN = 162


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--flank5", default="CAATCCGCCCTCACTACAACCGGGTCTCAAAGC")
    ap.add_argument("--flank3", default="TTCGAGAGACCCTACTCTGGCGTCGATGAGGGA")
    ap.add_argument("--intermediate-len", type=int, default=DEFAULT_INT_LEN)
    ap.add_argument("--strict", action="store_true",
                    help="exit non-zero if any check fails (default: warn and continue)")
    a = ap.parse_args()

    csv.field_size_limit(10 ** 8)
    rows = list(csv.DictReader(open(a.reference)))
    if not rows:
        sys.exit("reference is empty")
    need = {"gene", "order", "well", "vh_len_nt", "construct_dna"}
    missing = need - set(rows[0])
    if missing:
        sys.exit(f"reference is missing required columns: {sorted(missing)}")

    problems, meta = [], {}
    H, L = {}, {}
    seen = set()
    for r in rows:
        g, s = r["gene"], r["construct_dna"].strip().upper()
        if g in seen:
            problems.append(f"{g}: duplicate gene id")
            continue
        seen.add(g)
        if set(s) - set("ACGT"):
            problems.append(f"{g}: construct contains non-ACGT characters")
        if not s.startswith(a.flank5) or not s.endswith(a.flank3):
            problems.append(f"{g}: flanks do not match the expected 5'/3' flanks")
        if len(s) % 3 != 0:
            problems.append(f"{g}: length {len(s)} is not a multiple of 3")
        if r.get("length_nt") and int(r["length_nt"]) != len(s):
            problems.append(f"{g}: length_nt {r['length_nt']} != actual {len(s)}")
        vh = int(r["vh_len_nt"])
        il = int(r.get("intermediate_len_nt") or a.intermediate_len)
        cut = len(a.flank5) + vh + il // 2
        if not (0 < cut < len(s)):
            problems.append(f"{g}: half-split position {cut} is outside the construct")
            continue
        h, l = s[:cut], s[cut:]
        if h in H:
            problems.append(f"{g}: VH-side half is identical to {H[h]}")
        if l in L:
            problems.append(f"{g}: VL-side half is identical to {L[l]}")
        H[h], L[l] = g, g
        meta[g] = {"order": r["order"], "well": r["well"], "pool": r.get("pool", ""),
                   "framework": r.get("framework", ""), "vh_len": vh, "full_len": len(s),
                   "h_len": len(h), "l_len": len(l), "seq": s}

    with gzip.open("halves.fasta.gz", "wt") as fh:
        for g, m in meta.items():
            fh.write(f">H_{g}\n{m['seq'][:m['h_len']]}\n>L_{g}\n{m['seq'][m['h_len']:]}\n")
    with gzip.open("constructs.fasta.gz", "wt") as fh:
        for g, m in meta.items():
            fh.write(f">{g}\n{m['seq']}\n")
    for m in meta.values():
        m.pop("seq")
    json.dump(meta, open("gene_meta.json", "w"))

    wells = defaultdict(set)
    for m in meta.values():
        wells[m["order"]].add(m["well"])
    mates = Counter((m["order"], m["well"]) for m in meta.values())
    qc = {"n_constructs": len(meta), "n_problems": len(problems), "problems": problems[:200],
          "per_order": dict(Counter(m["order"] for m in meta.values())),
          "wells_per_order": {k: len(v) for k, v in wells.items()},
          "mates_per_well": dict(Counter(mates.values())),
          "frameworks": {o: sorted({m["framework"] for m in meta.values() if m["order"] == o})
                         for o in {m["order"] for m in meta.values()}},
          "length_nt": {"min": min(m["full_len"] for m in meta.values()),
                        "max": max(m["full_len"] for m in meta.values())},
          "half_split_rule": f"len(flank5) + vh_len_nt + intermediate_len_nt//2"}
    json.dump(qc, open("reference_qc.json", "w"), indent=1)
    print(json.dumps({k: qc[k] for k in ("n_constructs", "per_order", "wells_per_order",
                                         "mates_per_well", "frameworks", "length_nt", "n_problems")}, indent=1))
    for p in problems[:20]:
        print("PROBLEM:", p, file=sys.stderr)
    if problems and a.strict:
        sys.exit(f"{len(problems)} reference problems (--strict)")


if __name__ == "__main__":
    main()
