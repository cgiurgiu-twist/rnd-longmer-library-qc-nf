#!/usr/bin/env python3
"""
merge_metrics.py — fold the per-shard JSONs into the metric bundle every output is built from.

Denominators are named in the key, because mixing them is the easiest way to mislead:
  *_pct_of_all       share of every read in fastq_pass
  *_pct_of_assigned  share of reads with BOTH halves assigned (the pairing denominator)
  deletion rates     share of reads carrying BOTH primer sites and assigned to a variant

Writes metrics_core.json, metrics_extra.json, blockdel_summary.json, blockdel_scope.json,
per_gene.csv and per_well.csv.
"""
import argparse, csv, glob, json, os
from collections import Counter, defaultdict

try:
    import numpy as np
    _np_ok = True
except ImportError:          # numpy is not required for the merge itself
    _np_ok = False

PAIRED = ("clean_designed", "within_well", "cross_well", "cross_sublib")
ALL_CLS = PAIRED + ("one_half", "ambiguous", "unmapped")


def rates(c):
    tot = sum(v for k, v in c.items() if k in ALL_CLS)
    asg = sum(c.get(k, 0) for k in PAIRED)
    o = {"n_reads": tot, "n_both_halves_assigned": asg,
         "assigned_pct_of_all": round(100 * asg / tot, 4) if tot else None}
    for k in ALL_CLS:
        o[k + "_n"] = c.get(k, 0)
        o[k + "_pct_of_all"] = round(100 * c.get(k, 0) / tot, 4) if tot else None
    for k in PAIRED:
        o[k + "_pct_of_assigned"] = round(100 * c.get(k, 0) / asg, 4) if asg else None
    mis = sum(c.get(k, 0) for k in ("within_well", "cross_well", "cross_sublib"))
    o["total_mispaired_pct_of_assigned"] = round(100 * mis / asg, 4) if asg else None
    return o


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--classify", nargs="+", required=True)
    ap.add_argument("--blockdel", nargs="+", default=[])
    ap.add_argument("--meta", required=True)
    ap.add_argument("--outdir", default=".")
    a = ap.parse_args()
    meta = json.load(open(a.meta))
    os.makedirs(a.outdir, exist_ok=True)

    # ---------------- pairing ---------------------------------------------------------
    cls = Counter(); per_sub = defaultdict(Counter)
    per_gene = Counter(); per_well = defaultdict(Counter)
    rlen = Counter(); rlen_sub_cls = defaultdict(Counter); lenc = defaultdict(Counter)
    nfiles = 0
    for f in sorted(set(a.classify)):
        d = json.load(open(f)); nfiles += 1
        s = d["sublibrary"]
        for k, v in d["class_counts"].items():
            cls[k] += v; per_sub[s][k] += v
        for g, v in d["per_gene"].items():
            per_gene[g] += v
        for k, v in d["per_well"].items():
            for c, n in v.items():
                per_well[k][c] += n
        for b, v in d["rlen_hist_25bp"].items():
            rlen[int(b)] += v
        for c, h in d["rlen_hist_25bp_by_cls"].items():
            for b, v in h.items():
                rlen_sub_cls[f"{s}|{c}"][int(b)] += v
        for c, v in d["length_classes"].items():
            for k, n in v.items():
                lenc[c][k] += n

    total = sum(cls.get(k, 0) for k in ALL_CLS)
    core = {"total_reads": total, "n_files": nfiles,
            "overall": rates(cls), "per_sublib": {k: rates(v) for k, v in per_sub.items()},
            "integrity": {"F_full_length_reads": cls.get("_full_length_assigned", 0),
                          "F_pct_of_all": round(100 * cls.get("_full_length_assigned", 0) / total, 3) if total else None,
                          "P_of_F": cls.get("_full_length_clean", 0),
                          "P_over_F_pct": round(100 * cls.get("_full_length_clean", 0) / cls["_full_length_assigned"], 3)
                          if cls.get("_full_length_assigned") else None,
                          "note": "F = both halves assigned AND read >= 90% of designed length. "
                                  "(P and L)/F is an UPPER bound; L is measured separately."},
            "per_variant": {"designed": len(meta), "observed_ge1_read": len(per_gene),
                            "dropout_n": len(meta) - len(per_gene),
                            "dropout_pct": round(100 * (len(meta) - len(per_gene)) / len(meta), 5),
                            "note": "read the rarefaction/depth figures before quoting as a point estimate"},
            "n_wells_observed": len(per_well)}
    def binned_q(hist, qs=(0.05, 0.5, 0.95)):
        """Quantiles from the 25 bp read-length histogram. Binned, so exact to +/-12 nt —
        good enough for a length profile and it avoids keeping per-read lengths."""
        items = sorted(hist.items())
        tot = sum(v for _, v in items)
        if not tot:
            return [None] * len(qs)
        out, cum, i = [], 0, 0
        for q in qs:
            target = q * tot
            while i < len(items) and cum + items[i][1] < target:
                cum += items[i][1]; i += 1
            out.append(items[min(i, len(items) - 1)][0] + 12)
        return out

    rlen_cls = defaultdict(Counter)
    for k, v in rlen_sub_cls.items():
        c = k.split("|", 1)[1]
        for b, n in v.items():
            rlen_cls[c][b] += n
    extra_rl = {}
    for k, v in rlen_sub_cls.items():
        p5, med, p95 = binned_q(v)
        extra_rl[k] = [p5, med, p95]
    core["read_length_by_class"] = {}
    for c, h in rlen_cls.items():
        p5, med, p95 = binned_q(h)
        core["read_length_by_class"][c] = {"n": sum(h.values()), "p5": p5, "median": med, "p95": p95,
                                           "basis": "25 bp binned"}

    lclasses = {}
    for c, v in lenc.items():
        n = v["n"]
        lclasses[c] = {"n": n, "single_600_900": v["single_600_900"],
                       "pct_single": round(100 * v["single_600_900"] / n, 3) if n else None,
                       "fused_ge1200": v["fused_ge1200"],
                       "pct_fused": round(100 * v["fused_ge1200"] / n, 3) if n else None,
                       "mean_len": round(v["len_sum"] / n, 1) if n else None,
                       "median": core["read_length_by_class"].get(c, {}).get("median")}
    sm = Counter()
    for k, v in rlen_sub_cls.items():
        c = k.split("|", 1)[1]
        if c in PAIRED:
            sm[c] += sum(n for b, n in v.items() if 600 <= b <= 875)
    extra = {"rl_by_sublib_cls": extra_rl,
             "single_molecule_pairing": rates(sm),
             "single_molecule_pairing_note": "restricted to reads 600-900 nt, i.e. one molecule "
                                             "rather than a ligation-fused pair",
             "length_classes": lclasses,
             "rlen_hist_25bp": sorted([k, v] for k, v in rlen.items()),
             "rlen_hist_25bp_by_sublib_cls": {k: sorted([b, n] for b, n in v.items()) for k, v in rlen_sub_cls.items()}}
    # Rarefaction by multinomial thinning of the per-variant counts. For uniform random
    # subsampling of reads this is the exact distribution, so it needs no second pass over
    # the FASTQs; it is a model of subsampling, not a re-count, and is labelled as such.
    import random as _rnd
    counts = np.array([per_gene.get(g, 0) for g in meta], dtype=float) if _np_ok else None
    if _np_ok and counts.sum() > 0:
        p = counts / counts.sum()
        rar = []
        for frac in (0.05, 0.10, 0.25, 0.50, 0.75, 1.00):
            n = int(counts.sum() * frac)
            seen = int((1 - np.exp(-n * p)).sum().round()) if frac < 1 else int((counts > 0).sum())
            rar.append({"frac": frac, "n_reads": n, "n_variants": seen})
        core["rarefaction"] = rar
        core["rarefaction_basis"] = ("expected distinct variants under Poisson thinning of the "
                                     "observed per-variant counts; not a re-count of subsampled reads")
    json.dump(core, open(f"{a.outdir}/metrics_core.json", "w"), indent=1)
    json.dump(extra, open(f"{a.outdir}/metrics_extra.json", "w"), indent=1)
    with open(f"{a.outdir}/per_gene.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["gene", "n"])
        for g in sorted(meta):
            w.writerow([g, per_gene.get(g, 0)])
    with open(f"{a.outdir}/per_well.csv", "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["order", "well", "n_assigned", "n_clean", "n_within", "n_cross"])
        for k, v in sorted(per_well.items()):
            o, well = k.split("\t")
            w.writerow([o, well, sum(v.get(c, 0) for c in PAIRED), v.get("clean_designed", 0),
                        v.get("within_well", 0), v.get("cross_well", 0)])

    # ---------------- block deletion ---------------------------------------------------
    if a.blockdel:
        lad = Counter(); lad_sub = defaultdict(Counter); strand = defaultdict(Counter)
        elig = Counter(); elig_sub = defaultdict(Counter)
        calls = defaultdict(Counter); calls_sub = defaultdict(lambda: defaultdict(Counter))
        pg = defaultdict(lambda: [0, 0]); fl = Counter(); fl_sub = defaultdict(Counter)
        gs, gst, gen, ng = Counter(), Counter(), Counter(), Counter()
        idy = defaultdict(lambda: [0.0, 0]); conf = Counter(); nsub = 0; subsample = None
        ends = Counter()
        for f in sorted(set(a.blockdel)):
            d = json.load(open(f)); s = d["sublibrary"]; nsub += d["n_subsampled"]
            subsample = d["subsample"]
            for k, v in d["primer_ladder"].items():
                lad[k] += v; lad_sub[s][k] += v
            for k, v in d["strand"].items():
                for x, n in v.items():
                    strand[k][x] += n
            for k, v in d["eligible"].items():
                elig[k] += v; elig_sub[s][k] += v
            for thr, v in d["calls"].items():
                for k, n in v.items():
                    calls[thr][k] += n; calls_sub[s][thr][k] += n
            for g, (n, g50) in d["per_gene"].items():
                pg[g][0] += n; pg[g][1] += g50
            for k, v in d["full_length"].items():
                if k.startswith("cov_hist_") or k.startswith("covclean_hist_"):
                    fl[k] += v
                    if k.startswith("covclean_hist_"):
                        fl_sub[s][k] += v
                else:
                    fl[k] += v; fl_sub[s][k] += v
            for src, dst in ((d["gap_size_hist_10nt"], gs), (d["gap_start_hist_10nt"], gst), (d["gap_end_hist_10nt"], gen)):
                for b, v in src.items():
                    dst[int(b)] += v
            for k, v in d["n_gap_per_read"].items():
                ng[int(k)] += v
            for k in ("gap50", "no_gap"):
                if d["identity_median"][k] is not None:
                    idy[k][0] += d["identity_median"][k] * d["identity_n"][k]; idy[k][1] += d["identity_n"][k]
            conf["n"] += d["edlib_confirm"]["n"]; conf["ok"] += d["edlib_confirm"]["confirmed"]
            for k, v in d.get("endpoint_hists", {}).items():
                ends[k] += v
        ladtot = sum(lad.values()) or 1
        bd = {"n_subsampled": nsub, "subsample": subsample,
              "primer_ladder": dict(lad),
              "primer_ladder_pct": {k: round(100 * v / ladtot, 4) for k, v in lad.items()},
              "primer_ladder_by_sublib": {k: dict(v) for k, v in lad_sub.items()},
              "strand": {k: dict(v) for k, v in strand.items()},
              "eligible_n": elig["clean"],
              "block_deletion": {t: {"n": calls[t]["clean"], "pct": round(100 * calls[t]["clean"] / elig["clean"], 4)}
                                 for t in sorted(calls, key=int)} if elig["clean"] else {},
              "per_sublib": {s: {"n": elig_sub[s]["clean"], "g50": calls_sub[s]["50"]["clean"],
                                 "pct": round(100 * calls_sub[s]["50"]["clean"] / elig_sub[s]["clean"], 3)}
                             for s in elig_sub if elig_sub[s]["clean"]},
              "gap_size_hist_10nt": sorted([k, v] for k, v in gs.items()),
              "gap_start_hist_10nt": sorted([k, v] for k, v in gst.items()),
              "gap_end_hist_10nt": sorted([k, v] for k, v in gen.items()),
              "n_gap_per_read": {str(k): v for k, v in sorted(ng.items())},
              "identity": {k: (v[0] / v[1] if v[1] else None) for k, v in idy.items()},
              "edlib_confirm_pct": round(100 * conf["ok"] / conf["n"], 3) if conf["n"] else None}
        for what, key in (("r_st", "r_st_hist_10nt"), ("r_en", "r_en_hist_10nt")):
            h = defaultdict(list)
            for k, v in ends.items():
                kind, status, b = k.split("__")
                if kind == what:
                    h[status].append([int(b), v])
            bd[key] = {s: sorted(v) for s, v in h.items()}
        sc = {"n_subsampled": nsub, "subsample": subsample,
              "both_primer_n": lad.get("both", 0),
              "both_primer_pct": round(100 * lad.get("both", 0) / ladtot, 4),
              "eligible_any_cls_n": elig["all"],
              "blockdel50_n": calls["50"]["all"],
              "blockdel50_pct_of_eligible": round(100 * calls["50"]["all"] / elig["all"], 4) if elig["all"] else None,
              "per_sublib_any_cls": {s: {"n": elig_sub[s]["all"], "g50": calls_sub[s]["50"]["all"],
                                         "pct": round(100 * calls_sub[s]["50"]["all"] / elig_sub[s]["all"], 3)}
                                     for s in elig_sub if elig_sub[s]["all"]},
              "fl_matrix": {k: {d: (round(100 * fl[(("" if d == "both_primer" else d + "__")) + k] / fl["den_" + d], 3)
                                    if fl.get("den_" + d) else None)
                                for d in ("gene_assigned", "both_primer", "clean_designed")}
                            for k in ("aln_span_ge90", "readlen_ge90", "no_gap50", "readlen_and_nogap")},
              "fl_denominators": {d: fl.get("den_" + d, 0) for d in ("gene_assigned", "both_primer", "clean_designed")},
              "fl_readlen_nogap_by_sublib_both": {s: round(100 * fl_sub[s]["readlen_and_nogap"] / elig_sub[s]["all"], 3)
                                                  for s in elig_sub if elig_sub[s]["all"]},
              "cov_hist_10nt": sorted([int(k.split("_")[-1]), v] for k, v in fl.items() if k.startswith("cov_hist_")),
              "cov_hist_10nt_by_primer_status": {s: sorted([int(b), v] for kk, s2, b, v in
                                                           [(k.split("__")[0], k.split("__")[1], k.split("__")[2], v)
                                                            for k, v in ends.items() if k.startswith("cov__")]
                                                           if s2 == s)
                                                 for s in ("both", "fwd_only", "tail_only", "neither")},
              "cov_hist_10nt_clean_by_sublib": {s: sorted([int(k.split("_")[-1]), v] for k, v in fl_sub[s].items()
                                                          if k.startswith("covclean_hist_")) for s in fl_sub}}
        json.dump(bd, open(f"{a.outdir}/blockdel_summary.json", "w"), indent=1)
        json.dump(sc, open(f"{a.outdir}/blockdel_scope.json", "w"), indent=1)
        with open(f"{a.outdir}/blockdel_per_gene_all.csv", "w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["sublib", "gene", "n", "g50"])
            for g, (n, g50) in sorted(pg.items()):
                w.writerow([meta[g]["order"] if g in meta else "", g, n, g50])
        print(json.dumps({k: sc[k] for k in ("n_subsampled", "both_primer_pct", "eligible_any_cls_n",
                                             "blockdel50_pct_of_eligible", "per_sublib_any_cls")}, indent=1))
    print(json.dumps({"total_reads": core["total_reads"], "n_files": nfiles,
                      "overall_within_well_pct_of_assigned": core["overall"]["within_well_pct_of_assigned"],
                      "clean_pct_of_assigned": core["overall"]["clean_designed_pct_of_assigned"],
                      "dropout": core["per_variant"]["dropout_n"]}, indent=1))


if __name__ == "__main__":
    main()
