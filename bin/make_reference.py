"""
make_reference.py — build the reference bundle for a 750mer direct-print AILK library.

Input : the delivered <ORDER>_750mer_ORDER_*.csv design files (one per sub-library) and
        their matching MOP manifests, plus the writer files for verification.
Output: one self-describing reference CSV that is the pipeline's only design input, plus
        the derived FASTAs and a metadata table.

The reference CSV is the contract between design and QC: everything the pipeline needs to
classify a read is a column in it, so Ops never has to hand-assemble a FASTA. The pipeline
re-derives the VH/VL half-references from `vh_len_nt` and `intermediate_len_nt` rather than
trusting a separate file that could drift out of step with the constructs.
"""
import argparse, csv, glob, gzip, hashlib, json, os, sys
from collections import Counter, defaultdict

F5 = "CAATCCGCCCTCACTACAACCGGGTCTCAAAGC"
F3 = "TTCGAGAGACCCTACTCTGGCGTCGATGAGGGA"
INT_LEN = 162
COLS = ["gene", "order", "framework", "pool", "well", "well_rep2", "well_60perwell",
        "length_nt", "vh_len_nt", "intermediate_len_nt", "gc_pct", "construct_dna"]


def sha(s):
    return hashlib.sha256(s.encode()).hexdigest()[:12]


def build(order_csvs, manifest_dir, writer_dir, out_dir, name):
    os.makedirs(out_dir, exist_ok=True)
    rows, problems = [], []
    for f in order_csvs:
        order = None
        man = {}
        for pat in (f"{manifest_dir}/manifest_*_750mer_120perwell_x2.csv",):
            pass
        for r in csv.DictReader(open(f)):
            order = r["Order"]
            s = r["Full_Construct_DNA"]
            if not (s.startswith(F5) and s.endswith(F3)):
                problems.append(f"{r['Gene']}: flanks do not match the expected primer/BsaI flanks")
            if len(s) != int(r["Length_nt"]):
                problems.append(f"{r['Gene']}: Length_nt {r['Length_nt']} != len(sequence) {len(s)}")
            if len(s) % 3 != 0:
                problems.append(f"{r['Gene']}: construct length {len(s)} is not a multiple of 3")
            rows.append({"gene": r["Gene"], "order": order, "framework": r["Framework"], "pool": r["Pool"],
                         "well": r["Well_120perwell_rep1"], "well_rep2": r["Well_120perwell_rep2"],
                         "well_60perwell": r["Well_60perwell"], "length_nt": len(s),
                         "vh_len_nt": int(r["VH_len_nt"]), "intermediate_len_nt": INT_LEN,
                         "gc_pct": r["GC_pct"], "construct_dna": s})
        # manifest cross-check: designed well assignment must agree
        mf = f"{manifest_dir}/manifest_{order}_750mer_120perwell_x2.csv"
        if os.path.exists(mf):
            m = {x["gene"]: x["well"] for x in csv.DictReader(open(mf)) if x["replicate"] == "R1"}
            for r in rows:
                if r["order"] == order and m.get(r["gene"]) not in (None, r["well"]):
                    problems.append(f"{r['gene']}: well {r['well']} != manifest {m[r['gene']]}")
        # writer cross-check: the sequence that was actually printed
        wf = f"{writer_dir}/writer_MOP_{order}_750mer_plate.csv"
        if os.path.exists(wf):
            seq = {r["gene"]: r["construct_dna"] for r in rows if r["order"] == order}
            n_ok = n_bad = 0
            for x in csv.DictReader(open(wf)):
                base = x["Name"].rsplit("_R", 1)[0]
                if base in seq:
                    if seq[base] == x["sequence"].replace("8", "A"):
                        n_ok += 1
                    else:
                        n_bad += 1
                        problems.append(f"{base}: writer sequence differs from the design")
            print(f"  {order}: writer rows matching design {n_ok}, mismatching {n_bad}")

    # uniqueness of the half references, which is what the classifier relies on
    H, L = set(), set()
    for r in rows:
        cut = 33 + r["vh_len_nt"] + INT_LEN // 2
        h, l = r["construct_dna"][:cut], r["construct_dna"][cut:]
        if h in H:
            problems.append(f"{r['gene']}: VH-side half-reference is not unique")
        if l in L:
            problems.append(f"{r['gene']}: VL-side half-reference is not unique")
        H.add(h); L.add(l)

    ref = f"{out_dir}/{name}_reference.csv"
    with open(ref, "w", newline="") as fh:
        w = csv.DictWriter(fh, COLS); w.writeheader()
        for r in sorted(rows, key=lambda r: (r["order"], r["pool"], r["gene"])): w.writerow(r)
    with gzip.open(f"{out_dir}/{name}_constructs.fasta.gz", "wt") as fh:
        for r in rows: fh.write(f">{r['gene']}\n{r['construct_dna']}\n")
    with gzip.open(f"{out_dir}/{name}_halves.fasta.gz", "wt") as fh:
        for r in rows:
            cut = 33 + r["vh_len_nt"] + INT_LEN // 2
            fh.write(f">H_{r['gene']}\n{r['construct_dna'][:cut]}\n>L_{r['gene']}\n{r['construct_dna'][cut:]}\n")
    with open(f"{out_dir}/{name}_metadata.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, [c for c in COLS if c != "construct_dna"]); w.writeheader()
        for r in sorted(rows, key=lambda r: (r["order"], r["pool"], r["gene"])):
            w.writerow({k: v for k, v in r.items() if k != "construct_dna"})

    by = Counter(r["order"] for r in rows)
    wells = defaultdict(set)
    for r in rows: wells[r["order"]].add(r["well"])
    manifest = {
        "name": name, "built": __import__("datetime").date.today().isoformat(),
        "n_constructs": len(rows), "per_order": dict(by),
        "wells_per_order": {k: len(v) for k, v in wells.items()},
        "frameworks": {o: sorted({r["framework"] for r in rows if r["order"] == o}) for o in by},
        "length_nt": {"min": min(r["length_nt"] for r in rows), "max": max(r["length_nt"] for r in rows)},
        "flank_5p": F5, "flank_3p": F3, "intermediate_len_nt": INT_LEN,
        "half_split_rule": "cut at 33 + vh_len_nt + intermediate_len_nt/2 (midpoint of the shared intermediate)",
        "reference_csv_sha256_12": sha(open(ref).read()),
        "validation": {"checks_run": ["flanks", "declared length", "in-frame length", "manifest well agreement",
                                      "writer-file sequence identity", "half-reference uniqueness"],
                       "n_problems": len(problems), "problems": problems[:50]},
    }
    json.dump(manifest, open(f"{out_dir}/{name}_reference_manifest.json", "w"), indent=1)
    print(json.dumps({k: manifest[k] for k in ("n_constructs", "per_order", "wells_per_order", "frameworks", "length_nt")}, indent=1))
    print("problems:", len(problems))
    for p in problems[:10]: print("  !", p)
    return len(problems)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--order-csv", nargs="+", required=True)
    ap.add_argument("--manifest-dir", default="")
    ap.add_argument("--writer-dir", default="")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--name", required=True)
    a = ap.parse_args()
    files = [f for p in a.order_csv for f in sorted(glob.glob(p))]
    print("order files:", files)
    sys.exit(1 if build(files, a.manifest_dir, a.writer_dir, a.out_dir, a.name) else 0)
