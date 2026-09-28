#!/usr/bin/env python3
"""
run_controls.py — calibrate both callers on simulated reads before any rate is quoted.

No chimera or deletion rate from this pipeline is reportable without these numbers, so the
controls run on every execution against the library's OWN constructs rather than being
copied from a previous run — a reference with more well-mates, or shorter constructs, is a
harder problem for the classifier and the operating point has to be re-measured.

classifier
  negative : simulated clean reads at several error rates, uniform and clustered.
             Any chimera call here is a false positive; the reported ceiling is the 95%
             upper bound (3/n) when the count is zero.
  positive : injected within-well VH/VL swaps -> detection rate.

deletion caller
  truncation : reads with one end removed must be dropped by the both-primer filter.
  negative   : error-only reads, no programmed deletion -> false-call rate.
  positive   : programmed deletions of known size and position -> sensitivity and
               placement accuracy.
"""
import argparse, gzip, json, os, random, subprocess, sys, tempfile
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify_reads import classify_seq, label, load_aligner  # noqa: E402
from blockdel_reads import gaps_from_cigar, read_fq  # noqa: E402

COMP = str.maketrans("ACGT", "TGCA")


def rc(s):
    return s.translate(COMP)[::-1]


def mutate(seq, rate, rng):
    """ONT-like error: ~45% substitution, ~30% deletion, ~25% insertion."""
    out = []
    for c in seq:
        r = rng.random()
        if r >= rate:
            out.append(c)
        elif r < rate * .45:
            out.append(rng.choice("ACGT"))
        elif r < rate * .75:
            pass
        else:
            out.append(c); out.append(rng.choice("ACGT"))
    return "".join(out)


def clustered_mutate(seq, rate, rng, burst=12):
    """Same expected rate delivered in bursts. Nanopore error is not uniform, and a caller
    tuned only against uniform error will fabricate structural events in a burst."""
    s = list(seq)
    n_err = int(len(seq) * rate)
    done = 0
    while done < n_err:
        st = rng.randrange(max(1, len(s) - burst))
        for k in range(st, min(len(s), st + burst)):
            if done >= n_err:
                break
            r = rng.random()
            if r < .45:
                s[k] = rng.choice("ACGT")
            elif r < .75:
                s[k] = ""
            else:
                s[k] = s[k] + rng.choice("ACGT")
            done += 1
    return "".join(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--halves", required=True)
    ap.add_argument("--constructs", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n-negative", type=int, default=50000)
    ap.add_argument("--n-positive", type=int, default=500)
    ap.add_argument("--n-del", type=int, default=600)
    ap.add_argument("--seed", type=int, default=20260927)
    ap.add_argument("--flank5", default="CAATCCGCCCTCACTACAACCGGGTCTCAAAGC")
    ap.add_argument("--flank3", default="TTCGAGAGACCCTACTCTGGCGTCGATGAGGGA")
    a = ap.parse_args()
    import mappy

    rng = random.Random(a.seed)
    meta = json.load(open(a.meta))
    seqs = {}
    with gzip.open(a.constructs, "rt") as fh:
        name = None
        for line in fh:
            if line[0] == ">":
                name = line[1:].strip()
            else:
                seqs[name] = line.strip()
    aln = load_aligner(a.halves)
    mates = defaultdict(list)
    for g, m in meta.items():
        mates[(m["order"], m["well"])].append(g)
    genes = sorted(meta)
    out = {"classifier": {"negative": {}, "positive": {}}, "deletion_caller": {}}

    # ---------------- classifier negatives -------------------------------------------
    for cond, fn in (("flat", mutate), ("clustered", clustered_mutate)):
        for err in (0.01, 0.02, 0.03, 0.05, 0.08):
            n = a.n_negative if err <= 0.03 else max(2000, a.n_negative // 5)
            hit = 0
            for _ in range(n):
                g = rng.choice(genes)
                s = fn(seqs[g], err, rng)
                if rng.random() < .5:
                    s = rc(s)
                vh, vl, _ = classify_seq(aln, s)
                if label(vh, vl, meta) in ("within_well", "cross_well", "cross_sublib"):
                    hit += 1
            out["classifier"]["negative"][f"{cond}_{err:g}"] = {
                "n": n, "false_chimera": hit, "pct": round(100 * hit / n, 4),
                "upper95_pct_if_zero": round(100 * 3.0 / n, 4) if hit == 0 else None}
            print("neg", cond, err, out["classifier"]["negative"][f"{cond}_{err:g}"])

    # ---------------- classifier positives -------------------------------------------
    for err in (0.03, 0.05):
        hit = n = 0
        for _ in range(a.n_positive):
            gA = rng.choice(genes)
            m = meta[gA]
            pool = [x for x in mates[(m["order"], m["well"])] if x != gA]
            if not pool:
                continue
            gB = rng.choice(pool)
            s = mutate(seqs[gA][:m["h_len"]] + seqs[gB][meta[gB]["h_len"]:], err, rng)
            vh, vl, _ = classify_seq(aln, s)
            n += 1
            hit += label(vh, vl, meta) == "within_well"
        out["classifier"]["positive"][f"flat_{err:g}"] = {
            "n": n, "detected": hit, "pct": round(100 * hit / max(n, 1), 3)}
        print("pos", err, out["classifier"]["positive"][f"flat_{err:g}"])

    # ---------------- deletion caller ------------------------------------------------
    FWD, TAIL = a.flank5[:22], a.flank3[-22:]
    tmp = tempfile.mkdtemp()

    def through_filter(recs):
        src, dst, rest = f"{tmp}/i.fq", f"{tmp}/o.fq", f"{tmp}/r.fq"
        with open(src, "w") as fh:
            for i, s in enumerate(recs):
                fh.write(f"@r{i}\n{s}\n+\n{'I' * len(s)}\n")
        subprocess.run(["cutadapt", "-e", "0.15", "-O", "18", "--action=retain", "--quiet",
                        "--revcomp", "--untrimmed-output", rest, "-o", dst,
                        "-a", f"{FWD};required...{TAIL};required", src], check=True)
        return [s for _, s in read_fq(dst)]

    trunc = [mutate(seqs[(g := rng.choice(genes))][rng.randrange(60, 200):], 0.03, rng) for _ in range(500)]
    intact = [mutate(seqs[rng.choice(genes)], 0.03, rng) for _ in range(500)]
    out["deletion_caller"]["truncation_filter"] = {
        "truncated_in": len(trunc), "truncated_kept": len(through_filter(trunc)),
        "intact_in": len(intact), "intact_kept": len(through_filter(intact))}
    print("truncation", out["deletion_caller"]["truncation_filter"])

    def call_del(seq, gene, min_del=20):
        al = mappy.Aligner(seq=seqs[gene], preset="splice")
        best = None
        for h in al.map(seq):
            if best is None or h.mlen > best.mlen:
                best = h
        if best is None:
            return None
        gp = gaps_from_cigar(best, min_del)
        return max(gp, key=lambda x: x[0]) if gp else (0, 0)

    out["deletion_caller"]["negative"] = {}
    for cond, fn in (("flat", mutate), ("clustered", clustered_mutate)):
        for err in (0.01, 0.03, 0.05, 0.08):
            n = false = 0
            for _ in range(a.n_del):
                g = rng.choice(genes)
                s = fn(seqs[g], err, rng)
                kept = through_filter([s])
                if not kept:
                    continue
                n += 1
                c = call_del(kept[0], g)
                if c and c[0] >= 50:
                    false += 1
            out["deletion_caller"]["negative"][f"{cond}_{err:g}"] = {
                "n_after_filter": n, "false_calls": false, "pct": round(100 * false / max(n, 1), 4),
                "upper95_pct_if_zero": round(100 * 3.0 / max(n, 1), 4) if false == 0 else None}
            print("del neg", cond, err, out["deletion_caller"]["negative"][f"{cond}_{err:g}"])

    out["deletion_caller"]["positive"] = {}
    for size in (20, 30, 50, 75, 100, 150, 200, 300):
        n = det = 0
        poserr = []
        for _ in range(a.n_del):
            g = rng.choice(genes)
            s0 = seqs[g]
            lo, hi = 40, len(s0) - 40 - size
            if hi <= lo:
                continue
            st = rng.randrange(lo, hi)
            s = mutate(s0[:st] + s0[st + size:], 0.03, rng)
            kept = through_filter([s])
            if not kept:
                continue
            n += 1
            c = call_del(kept[0], g)
            if c and c[0] >= min(size, 50) - 5:
                det += 1
                poserr.append(abs(c[1] - st))
        poserr.sort()
        out["deletion_caller"]["positive"][str(size)] = {
            "n_after_filter": n, "detected": det, "sensitivity_pct": round(100 * det / max(n, 1), 2),
            "median_position_error_nt": poserr[len(poserr) // 2] if poserr else None}
        print("del pos", size, out["deletion_caller"]["positive"][str(size)])

    json.dump(out, open(a.out, "w"), indent=1)


if __name__ == "__main__":
    main()
