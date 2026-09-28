#!/usr/bin/env python3
"""
analyze.py — derived statistics for a 750mer direct-print library QC run.

Turns the merged metric bundle into derived.json: per-sub-library uniformity, per-well
depth and chimera, the block-deletion anatomy and junction analysis, the primer-collision
attribution, and per-sub-library screening depth.

Every number in the report, workbook, Slack summary and customer collateral is read from
this file or from the merged metric JSONs — nothing downstream recomputes a rate, so a
figure and a sentence can never disagree.

Bounds carried through deliberately: mosaic is a LOWER bound, intact fractions are UPPER
bounds, and screening-depth cell counts are MINIMA.
"""
import csv, json, math, random, collections, os
import numpy as np

import argparse
AP = argparse.ArgumentParser(description="Derived QC statistics from the metric bundle.")
AP.add_argument("--metrics-dir", default=".")
AP.add_argument("--reference", required=True, help="design reference CSV (for construct sequences)")
AP.add_argument("--meta", required=True, help="gene_meta.json from build_reference")
AP.add_argument("--microhomology", default=None, help="optional per_well_microhomology.csv")
AP.add_argument("--design-screen", default=None, help="optional design_primer_screen.json")
AP.add_argument("--out", default=None)
AP.add_argument("--seed", type=int, default=20260927)
A = AP.parse_args()
D = A.metrics_dir.rstrip("/") + "/"
RNG = random.Random(A.seed)
J = lambda n: json.load(open(D + n))
csv.field_size_limit(10 ** 8)
_META = json.load(open(A.meta))
ORD = sorted({m["order"] for m in _META.values()})


def meta():
    return _META


def constructs():
    out = {}
    for r in csv.DictReader(open(A.reference)):
        out[r["gene"]] = (r["construct_dna"], int(r["vh_len_nt"]))
    return out


def uniformity(counts, designed):
    n = np.array(counts, float)
    n = np.concatenate([n, np.zeros(designed - len(n))])
    med = float(np.median(n))
    srt = np.sort(n); k = np.arange(1, len(srt) + 1)
    gini = float(2 * np.sum(k * srt) / (len(srt) * srt.sum()) - (len(srt) + 1) / len(srt))
    return dict(n_designed=designed, n_observed=int((n > 0).sum()), dropout=int((n == 0).sum()),
                dropout_pct=round(100 * float((n == 0).mean()), 4),
                min=float(n.min()), p5=float(np.percentile(n, 5)), p25=float(np.percentile(n, 25)),
                median=med, p75=float(np.percentile(n, 75)), p95=float(np.percentile(n, 95)),
                mean=float(n.mean()), cv=float(n.std() / n.mean()), gini=gini,
                p95p5=(float(np.percentile(n, 95) / np.percentile(n, 5)) if np.percentile(n, 5) > 0 else None),
                p95p5_note=(None if np.percentile(n, 5) > 0 else
                            "undefined: the 5th percentile is 0, i.e. >=5% of variants have no reads"),
                within2x=float(100 * np.mean((n >= 0.5 * med) & (n <= 2 * med))),
                n_below_10=int((n < 10).sum()))


CELL_CAP = 1e10


def solve(target, p, f):
    """Cells needed for `target`% of variants to be seen at least once, Poisson.

    Returns None rather than the search ceiling when the target is unreachable — a variant
    with zero abundance or zero intact fraction can never be recovered, and reporting the
    bisection bound as a cell count would be a fabricated number."""
    f = np.asarray(f, dtype=float)
    eff = np.asarray(p, dtype=float) * f          # per-variant rate of an intact, correct clone
    if not np.all(np.isfinite(eff)) or eff.max() <= 0:
        return None
    if 100 * (eff <= 0).mean() > (100 - target):
        return None                                # target unreachable: too many dead variants
    lo, hi = 1e2, CELL_CAP
    for _ in range(200):
        mid = (lo + hi) / 2
        if np.mean(1 - np.exp(-mid * p * f)) * 100 < target:
            lo = mid
        else:
            hi = mid
    v = (lo + hi) / 2
    return None if v >= CELL_CAP * 0.99 else v


def repeat_len(seq, a, b, cap=40):
    d = 0
    while d < cap and a + d < len(seq) and b + d < len(seq) and seq[a + d] == seq[b + d]: d += 1
    u = 0
    while u < cap and a - 1 - u >= 0 and b - 1 - u >= 0 and seq[a - 1 - u] == seq[b - 1 - u]: u += 1
    return u + d


def pearson_ci(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    r = float(np.corrcoef(x, y)[0, 1]); n = len(x)
    z = 0.5 * math.log((1 + r) / (1 - r)); se = 1 / math.sqrt(n - 3)
    return dict(r=round(r, 4), ci=[round(math.tanh(z - 1.96 * se), 4), round(math.tanh(z + 1.96 * se), 4)], n=n)


def main():
    M = meta(); cons = constructs()
    core = J("metrics_core.json"); ext = J("metrics_extra.json")
    BD = J("blockdel_summary.json"); SC = J("blockdel_scope.json")
    out = {}

    # ---------------- uniformity, per sub-library and pooled --------------------------
    pg = {r["gene"]: int(r["n"]) for r in csv.DictReader(open(D + "per_gene.csv"))}
    uni = {}
    for o in ORD:
        genes = [g for g in M if M[g]["order"] == o]
        uni[o] = uniformity([pg[g] for g in genes if g in pg], len(genes))
        uni[o]["dropped_genes"] = sorted(g for g in genes if g not in pg)[:50]
    uni["ALL"] = uniformity(list(pg.values()), len(M))
    out["uniformity"] = uni

    # per-well depth and chimera ----------------------------------------------------
    pw = list(csv.DictReader(open(D + "per_well.csv")))
    wells = {}
    for r in pw:
        k = (r["order"], r["well"])
        wells[k] = dict(n_assigned=int(r["n_assigned"]), n_clean=int(r["n_clean"]),
                        n_within=int(r["n_within"]), n_cross=int(r["n_cross"]))
    pwstat = {}
    for o in ORD:
        ks = [k for k in wells if k[0] == o]
        wi = np.array([100 * wells[k]["n_within"] / wells[k]["n_assigned"] for k in ks])
        cr = np.array([100 * wells[k]["n_cross"] / wells[k]["n_assigned"] for k in ks])
        dep = np.array([wells[k]["n_clean"] for k in ks], float)
        nd = len({M[g]["well"] for g in M if M[g]["order"] == o})
        pwstat[o] = dict(n_wells_designed=nd, n_wells_observed=len(ks),
                         within_median=float(np.median(wi)), within_max=float(wi.max()),
                         cross_median=float(np.median(cr)), cross_max=float(cr.max()),
                         depth_median=float(np.median(dep)), depth_min=float(dep.min()),
                         depth_max=float(dep.max()), depth_cv=float(dep.std() / dep.mean()))
    out["per_well"] = pwstat

    # microhomology regression — WITHIN sub-library (AILK013 lesson: pooled r is Simpson-prone)
    mh = {(r["order"], r["well"]): r for r in csv.DictReader(open(A.microhomology))} if A.microhomology else {}
    reg = {}
    for scope in ORD + ["pooled"]:
        xs, ys = collections.defaultdict(list), []
        for k, w in wells.items():
            if scope != "pooled" and k[0] != scope: continue
            m = mh.get(k)
            if not m or int(m["n_mates"]) != 120: continue
            ys.append(100 * w["n_within"] / w["n_assigned"])
            for f in ("shared_18mers", "shared_25mers", "max_lcs", "n_pairs_ge_40"):
                xs[f].append(float(m[f]))
        reg[scope] = {f: pearson_ci(v, ys) for f, v in xs.items() if len(v) > 5}
        reg[scope]["mean_within_events_per_well"] = float(np.mean(
            [w["n_within"] for k, w in wells.items() if scope == "pooled" or k[0] == scope]))
    out["regression"] = reg
    if mh:
        out["microhomology_design"] = {o: {f: float(np.median([float(r[f]) for k, r in mh.items() if k[0] == o and r.get(f)]))
                                           for f in ("shared_18mers", "shared_25mers", "max_lcs", "n_pairs_ge_40")} for o in ORD}

    # ---------------- block deletions ---------------------------------------------------
    ev = [e for e in csv.DictReader(open(D + "blockdel_events.csv")) if e["gene"] in cons]
    sz = np.array([int(e["gap_len"]) for e in ev], float)
    bd = {"n_events_sampled": len(ev), "size_median": float(np.median(sz)),
          "size_iqr": [float(np.percentile(sz, 25)), float(np.percentile(sz, 75))],
          "one_gap_pct": 100 * BD["n_gap_per_read"].get("1", 0) / sum(BD["n_gap_per_read"].values())}
    bd["size_by_order"] = {}
    for o in ORD:
        s = np.array([int(e["gap_len"]) for e in ev if M[e["gene"]]["order"] == o], float)
        bd["size_by_order"][o] = dict(n=len(s), median=float(np.median(s)) if len(s) else None)
    # anatomy
    seg = collections.Counter()
    for e in ev:
        seq, vh = cons[e["gene"]]; st, ln = int(e["gap_start"]), int(e["gap_len"])
        b = [33, 33 + vh, 33 + vh + 162, len(seq) - 33, len(seq)]  # flank5 | VH | intermediate | VL | flank3
        w = lambda p: next(i for i, x in enumerate(b) if p < x) if p < b[-1] else 4
        a, z = w(st), w(st + ln - 1)
        seg["spans" if a != z else ["5' flank", "VH", "intermediate", "VL", "3' flank"][a]] += 1
    tot = sum(seg.values())
    bd["anatomy"] = {k: {"n": v, "pct": round(100 * v / tot, 2)} for k, v in seg.items()}
    # microhomology at the junction, both directions, vs position-matched null
    obs, null = [], []
    for e in ev:
        s = cons[e["gene"]][0]; a, ln = int(e["gap_start"]), int(e["gap_len"])
        if a + ln >= len(s): continue
        obs.append(repeat_len(s, a, a + ln))
        lo, hi = 33, len(s) - 33 - ln
        if hi > lo:
            a2 = RNG.randint(lo, hi); null.append(repeat_len(s, a2, a2 + ln))
    obs, null = np.array(obs), np.array(null)
    bd["junction_repeat"] = dict(obs_mean=float(obs.mean()), null_mean=float(null.mean()),
                                 obs_pct_ge10=float(100 * (obs >= 10).mean()),
                                 null_pct_ge10=float(100 * (null >= 10).mean()))
    # reproducibility of junction per construct
    by = collections.defaultdict(list)
    for e in ev: by[e["gene"]].append(int(e["gap_start"]))
    modal, nmodal = [], []
    for g, v in by.items():
        if len(v) < 10: continue
        c = collections.Counter(x // 10 for x in v); modal.append(c.most_common(1)[0][1] / len(v))
        c2 = collections.Counter(RNG.randrange(60) for _ in v); nmodal.append(c2.most_common(1)[0][1] / len(v))
    bd["reproducibility"] = dict(n_genes=len(modal), modal_median=float(np.median(modal)) if modal else None,
                                 null_median=float(np.median(nmodal)) if nmodal else None)
    # primer-collision product: junction end at 175-195 AND start within the 5' flank region (<40)
    att = {}
    for o in ORD:
        e = [x for x in ev if x["gene"] in M and M[x["gene"]]["order"] == o]
        if not e:
            att[o] = {"n_events": 0, "at_primer_site": 0, "pct": None, "pct_points": None,
                      "pct_points_clean": None, "note": "no sampled deletion events in this sub-library"}
            continue
        at = [x for x in e if 170 <= int(x["gap_start"]) + int(x["gap_len"]) <= 200 and int(x["gap_start"]) <= 40]
        sizes = [int(x["gap_len"]) for x in at]
        att[o] = dict(n_events=len(e), at_primer_site=len(at), pct=round(100 * len(at) / max(1, len(e)), 2),
                      median_size=float(np.median(sizes)) if sizes else None,
                      start_mode=collections.Counter(int(x["gap_start"]) for x in at).most_common(3),
                      end_mode=collections.Counter(int(x["gap_start"]) + int(x["gap_len"]) for x in at).most_common(3))
        rate = SC["per_sublib_any_cls"].get(o, {}).get("pct", 0.0)
        att[o]["pct_points"] = round(rate * att[o]["pct"] / 100, 2)
        att[o]["pct_points_clean"] = round(BD["per_sublib"].get(o, {}).get("pct", 0.0) * att[o]["pct"] / 100, 2)
    bd["primer_site_attribution"] = att
    # the non-primer-site population: start/end medians per order
    rest = {}
    for o in ORD:
        e = [x for x in ev if x["gene"] in M and M[x["gene"]]["order"] == o
             and not (170 <= int(x["gap_start"]) + int(x["gap_len"]) <= 200 and int(x["gap_start"]) <= 40)]
        rest[o] = dict(n=len(e),
                       start_median=float(np.median([int(x["gap_start"]) for x in e])) if e else None,
                       end_median=float(np.median([int(x["gap_start"]) + int(x["gap_len"]) for x in e])) if e else None,
                       size_median=float(np.median([int(x["gap_len"]) for x in e])) if e else None)
    bd["non_primer_population"] = rest
    # joint start/end 20x20 bins: largest bins per order
    top = {}
    for o in ORD:
        c = collections.Counter((int(x["gap_start"]) // 20 * 20, (int(x["gap_start"]) + int(x["gap_len"])) // 20 * 20)
                                for x in ev if x["gene"] in M and M[x["gene"]]["order"] == o)
        top[o] = [[k[0], k[1], v] for k, v in c.most_common(4)]
    bd["top_joint_bins"] = top

    # per-gene and per-well deletion spread, per order (both-primer, gene-assigned reads)
    pga = list(csv.DictReader(open(D + "blockdel_per_gene_all.csv")))
    spread = {}; pwdel = {}
    for o in ORD:
        rows = [r for r in pga if r["sublib"] == o and r["gene"] in M and M[r["gene"]]["order"] == o]
        big = [r for r in rows if int(r["n"]) >= 50]
        rate = np.array([int(r["g50"]) / int(r["n"]) for r in big])
        spread[o] = (dict(n_genes=len(big), p5=float(np.percentile(rate, 5) * 100), median=float(np.median(rate) * 100),
                          p95=float(np.percentile(rate, 95) * 100), pct_under5=float(100 * (rate < .05).mean()),
                          pct_over50=float(100 * (rate > .5).mean()))
                     if len(big) else dict(n_genes=0, note="no variant reached the 50-read minimum in this sub-library"))
        w = collections.defaultdict(lambda: [0, 0])
        for r in rows:
            k = M[r["gene"]]["well"]; w[k][0] += int(r["n"]); w[k][1] += int(r["g50"])
        wr = np.array([100 * b / a for a, b in w.values() if a >= 200])
        pwdel[o] = (dict(n_wells=len(wr), min=float(wr.min()), median=float(np.median(wr)), max=float(wr.max()))
                    if len(wr) else dict(n_wells=0, note="no well reached the 200-read minimum in this sub-library"))
        pwdel[o]["_wells"] = {k: round(100 * b / a, 3) for k, (a, b) in w.items() if a >= 200}
    bd["per_gene_spread"] = spread
    bd["per_well"] = pwdel
    # reads from gene-level g50 in sublib that do not belong to it (sublib mislabel check)
    bd["foreign_gene_rows"] = sum(int(r["n"]) for r in pga if r["gene"] in M and M[r["gene"]]["order"] != r["sublib"])
    out["blockdel"] = bd

    # ---------------- recoverability, per sub-library ----------------------------------
    rec = {}
    ga = {(r["sublib"], r["gene"]): (int(r["n"]), int(r["g50"])) for r in pga}
    for o in ORD:
        genes = [g for g in M if M[g]["order"] == o]
        cnt = np.array([pg.get(g, 0) for g in genes], float)
        p = cnt / cnt.sum()
        P = core["per_sublib"][o]["clean_designed_pct_of_assigned"] / 100
        fr = []
        orate = SC["per_sublib_any_cls"].get(o, {}).get("pct", 0.0) / 100
        for g in genes:
            n, d = ga.get((o, g), (0, 0))
            fr.append(1 - (d / n if n >= 20 else orate))
        fr = np.array(fr)
        r = {"P": P, "mass_weighted_intact": float(np.sum(p * fr)), "variants_below_50pct_intact": int((fr < .5).sum())}
        for t in (90, 95, 99):
            a_, b_ = solve(t, p, P), solve(t, p, P * fr)
            r[str(t)] = dict(ignoring_deletions=int(round(a_)) if a_ else None,
                             counting_deletions=int(round(b_)) if b_ else None,
                             ratio=round(b_ / a_, 3) if (a_ and b_) else None,
                             note=None if (a_ and b_) else
                             "unreachable: too many variants have no reads at this depth")
        rec[o] = r
    out["recoverability"] = rec

    # ---------------- design side ----------------------------------------------------------
    if A.design_screen:
        out["design"] = json.load(open(A.design_screen))
    json.dump(out, open(A.out or (D + "derived.json"), "w"), indent=1, default=float)
    print(json.dumps({k: out[k] for k in ("uniformity",)}, default=float)[:1500])
    print(json.dumps(out["blockdel"]["primer_site_attribution"]))
    print(json.dumps(out["recoverability"]))


if __name__ == "__main__":
    main()
