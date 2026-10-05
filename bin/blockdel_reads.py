#!/usr/bin/env python3
"""
blockdel_reads.py — internal block deletions, primer ladder and strand, for one FASTQ shard.

A block deletion is a molecule that still carries BOTH ends — both primer binding sites —
but is missing a stretch from the middle. The read has to prove it still has both ends
before anything can be said about its middle, so:

  1. cutadapt with a linked adapter requiring both primers (--action=retain, so the read
     still spans the design). Reads that fail are end-truncated; they are profiled
     separately rather than counted as deletions. A second trim of the same linked pair
     on the both-primer reads gives the primer-to-primer molecule length used for the
     full-length 90% check (raw FASTQ length still includes ONT adapter/end-prep).
  2. minimap2 splice preset against the read's OWN assigned construct, never against an
     index of every construct — a pilot that took the best hit across all constructs
     called a gap in 13% of reads at 0.75 identity, because the splice model bridges a
     near-neighbour's CDR mismatches as false introns.
  3. a call is a single N/D operation >= --call-del nt, cross-checked with an edlib global
     alignment that has no splice model at all.

Emits one compact JSON per shard plus a capped sample of individual events (for the
breakpoint-position and junction analyses downstream).
"""
import argparse, gzip, json, os, shutil, subprocess, sys, tempfile
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from classify_reads import classify_seq, label, load_aligner, open_fastq  # noqa: E402


def read_fq(path):
    out = []
    with open(path) as fh:
        while True:
            h = fh.readline()
            if not h:
                break
            s = fh.readline().strip(); fh.readline(); fh.readline()
            out.append((h[1:].rstrip("\n"), s))
    return out


def gaps_from_cigar(hit, min_del):
    g, rpos = [], hit.r_st
    for ln, op in hit.cigar:
        if op in (0, 7, 8):
            rpos += ln
        elif op in (2, 3):
            if ln >= min_del:
                g.append((ln, rpos))
            rpos += ln
    return g


def edlib_gaps(edlib, ref, seq, min_del, bridge=12):
    """Splice-model-free confirmation. edlib has no gap-open cost, so one long deletion is
    split by chance identities inside it; operations separated by < bridge matches merge."""
    cig = edlib.align(seq, ref, mode="NW", task="path")["cigar"]
    raw, rpos, i = [], 0, 0
    while i < len(cig):
        j = i
        while cig[j].isdigit():
            j += 1
        ln, op = int(cig[i:j]), cig[j]
        if op in "=XM":
            rpos += ln
        elif op == "D":
            raw.append([rpos, ln]); rpos += ln
        i = j + 1
    m = []
    for st, ln in raw:
        if m and st - (m[-1][0] + m[-1][1]) <= bridge:
            m[-1][1] = st + ln - m[-1][0]
        else:
            m.append([st, ln])
    return [(ln, st) for st, ln in m if ln >= min_del]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--halves", required=True)
    ap.add_argument("--constructs", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--sublibrary", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--events-out", required=True)
    ap.add_argument("--subsample", type=int, default=8, help="keep every Nth read")
    ap.add_argument("--min-del", type=int, default=20, help="record gaps from this size")
    ap.add_argument("--call-del", type=int, default=50, help="report gaps from this size")
    ap.add_argument("--max-events", type=int, default=100000)
    ap.add_argument("--flank5", default="CAATCCGCCCTCACTACAACCGGGTCTCAAAGC")
    ap.add_argument("--flank3", default="TTCGAGAGACCCTACTCTGGCGTCGATGAGGGA")
    ap.add_argument("--cpus", type=int, default=1)
    a = ap.parse_args()
    import mappy, edlib

    FWD, TAIL = a.flank5[:22], a.flank3[-22:]
    meta = json.load(open(a.meta))
    refs = {}
    with gzip.open(a.constructs, "rt") as fh:
        name = None
        for line in fh:
            if line[0] == ">":
                name = line[1:].strip()
            else:
                refs[name] = line.strip()
    aln = load_aligner(a.halves)
    tmp = tempfile.mkdtemp()
    raw = f"{tmp}/raw.fq"
    n = 0
    with open_fastq(a.fastq) as fh, open(raw, "w") as out:
        for i, line in enumerate(fh):
            if i % 4 != 1:
                continue
            if (i // 4) % a.subsample:
                continue
            s = line.strip()
            if len(s) < 60:
                continue
            out.write(f"@r{n}\n{s}\n+\n{'I' * len(s)}\n")
            n += 1

    def cutadapt(src, dst, args):
        subprocess.run(["cutadapt", "-e", "0.15", "-O", "18", "--action=retain", "--quiet",
                        "-j", str(a.cpus), "-o", dst] + args + [src], check=True)

    both, rest = f"{tmp}/both.fq", f"{tmp}/rest.fq"
    f_only, t_only = f"{tmp}/f.fq", f"{tmp}/t.fq"
    span = {}
    n_span_fail = 0
    if n:
        cutadapt(raw, both, ["--revcomp", "--untrimmed-output", rest,
                             "-a", f"{FWD};required...{TAIL};required"])
        cutadapt(rest, f_only, ["--revcomp", "--discard-untrimmed", "-g", FWD])
        cutadapt(rest, t_only, ["--revcomp", "--discard-untrimmed", "-a", TAIL])
        fset = {h.split()[0] for h, _ in read_fq(f_only)}
        tset = {h.split()[0] for h, _ in read_fq(t_only)}
        recs = [(h, s, "both") for h, s in read_fq(both)]
        for h, s in read_fq(rest):
            rid = h.split()[0]
            recs.append((h, s, "fwd_only" if rid in fset else ("tail_only" if rid in tset else "neither")))
        # Second pass: trim (do not retain) the already-oriented both-primer reads so the
        # leftover sequence is the interior between the two primer matches. Primer-span
        # length is interior + both primer matches; ONT adapter/end-prep outside the
        # primers is dropped. --revcomp is off because the retain pass already oriented.
        both_trim = f"{tmp}/both_trim.fq"
        subprocess.run(["cutadapt", "-e", "0.15", "-O", "18", "--quiet", "--discard-untrimmed",
                        "-j", str(a.cpus), "-o", both_trim, "-a", f"{FWD}...{TAIL}", both],
                       check=True)
        for h, s in read_fq(both_trim):
            span[h.split()[0]] = len(s) + len(FWD) + len(TAIL)
    else:
        recs = []

    ladder, strand, ends = Counter(), defaultdict(Counter), Counter()
    elig = Counter()                      # both-primer + gene-assigned, by sublib
    calls = defaultdict(Counter)          # threshold -> counts
    per_gene = defaultdict(lambda: [0, 0])
    gap_size, gap_start, gap_end = Counter(), Counter(), Counter()
    ngap = Counter()
    fl = Counter()
    mol_hist, mol_hist_gene, mol_delta, raw_both = Counter(), Counter(), Counter(), Counter()
    idy_gap, idy_nogap = [], []
    events, confirmed, confirm_tot = [], 0, 0
    cache = {}

    for h, s, status in recs:
        ladder[status] += 1
        rid = h.split()[0]
        mol = span.get(rid) if status == "both" else None
        if status == "both":
            raw_both[min(len(s), 1200) // 10 * 10] += 1
            if mol is None:
                n_span_fail += 1
            else:
                mol_hist[min(mol, 1200) // 10 * 10] += 1
        flipped = h.rstrip().endswith("rc")
        vh, vl, _ = classify_seq(aln, s)
        cl = label(vh, vl, meta)
        gene = vh if (vh and vl and vh == vl) else ((vh or vl) if not (vh and vl) else None)
        if status == "both":
            strand["both"]["reverse" if flipped else "forward"] += 1
        if gene is None or gene not in refs:
            continue
        ref = refs[gene]
        al = cache.get(gene)
        if al is None:
            if len(cache) > 512:
                cache.clear()
            al = cache[gene] = mappy.Aligner(seq=ref, preset="splice")
        best = None
        for hit in al.map(s):
            if best is None or hit.mlen > best.mlen:
                best = hit
        if best is None:
            continue
        if status != "both":
            strand[status]["forward" if best.strand == 1 else "reverse"] += 1
        # endpoint profiles for every primer class: a synthesis-truncation series shares one
        # end and varies the other; pore dropout varies both.
        ends[f"r_st__{status}__{min(best.r_st, 750) // 10 * 10}"] += 1
        ends[f"r_en__{status}__{min(best.r_en, 750) // 10 * 10}"] += 1
        gp = gaps_from_cigar(best, a.min_del)
        mx = max((x[0] for x in gp), default=0)
        gapsum = sum(x[0] for x in gp)
        aln_len = best.r_en - best.r_st
        cov = aln_len - gapsum
        full = meta[gene]["full_len"]
        ends[f"cov__{status}__{min(cov, 800) // 10 * 10}"] += 1

        # Full-length flags are counted against three denominators, because the three
        # reasonable definitions of "full length" disagree by >10 points on the same reads
        # and the denominator has to be stated with the number. Alignment span straddles an
        # internal deletion, so it cannot see this defect at all — it is reported to show that.
        # readlen_* uses primer-to-primer span (both primers required), not raw FASTQ length.
        readlen_ok = mol is not None and mol >= 0.9 * full
        flags = {"aln_span_ge90": aln_len >= 0.9 * full,
                 "readlen_ge90": readlen_ok,
                 "no_gap50": mx < a.call_del,
                 "readlen_and_nogap": readlen_ok and mx < a.call_del}
        if mol is not None:
            mol_hist_gene[min(mol, 1200) // 10 * 10] += 1
            dlt = int(mol - full)
            mol_delta[max(-200, min(200, dlt // 10 * 10))] += 1
        fl["den_gene_assigned"] += 1
        for k, ok in flags.items():
            if ok:
                fl["gene_assigned__" + k] += 1
        if cl == "clean_designed":
            fl["den_clean_designed"] += 1
            for k, ok in flags.items():
                if ok:
                    fl["clean_designed__" + k] += 1
        if status != "both":
            continue
        elig["all"] += 1
        if cl == "clean_designed":
            elig["clean"] += 1
        for thr in (20, 30, 50, 100, 200):
            if mx >= thr:
                calls[thr]["all"] += 1
                if cl == "clean_designed":
                    calls[thr]["clean"] += 1
        per_gene[gene][0] += 1
        if mx >= a.call_del:
            per_gene[gene][1] += 1
            ln, st = max(gp, key=lambda x: x[0])
            gap_size[min(ln, 600) // 10 * 10] += 1
            gap_start[min(st, 750) // 10 * 10] += 1
            gap_end[min(st + ln, 750) // 10 * 10] += 1
            ngap[len(gp)] += 1
            if aln_len - gapsum > 0:
                idy_gap.append(best.mlen / (aln_len - gapsum))
            confirm_tot += 1
            try:
                confirmed += 1 if edlib_gaps(edlib, ref, s, a.min_del) else 0
            except Exception:
                pass
            if len(events) < a.max_events:
                events.append((a.sublibrary, gene, st, ln, len(s), mol if mol is not None else "", len(ref)))
        elif aln_len - gapsum > 0:
            idy_nogap.append(best.mlen / (aln_len - gapsum))
        fl["den_both_primer"] += 1
        for k, ok in flags.items():
            if ok:
                fl[k] += 1
        fl["cov_hist_" + str(min(cov, 800) // 10 * 10)] += 1
        if cl == "clean_designed":
            fl["covclean_hist_" + str(min(cov, 800) // 10 * 10)] += 1

    def med(v):
        v = sorted(v)
        return v[len(v) // 2] if v else None

    out = {"sublibrary": a.sublibrary, "fastq": os.path.basename(a.fastq),
           "subsample": a.subsample, "n_subsampled": n,
           "primer_ladder": dict(ladder), "strand": {k: dict(v) for k, v in strand.items()},
           "endpoint_hists": dict(ends),
           "eligible": dict(elig), "calls": {str(k): dict(v) for k, v in calls.items()},
           "per_gene": {g: v for g, v in per_gene.items()},
           "gap_size_hist_10nt": dict(gap_size), "gap_start_hist_10nt": dict(gap_start),
           "gap_end_hist_10nt": dict(gap_end), "n_gap_per_read": dict(ngap),
           "full_length": dict(fl),
           "mol_len_hist_10nt": dict(mol_hist),
           "mol_len_hist_10nt_gene": dict(mol_hist_gene),
           "mol_minus_design_hist_10nt": dict(mol_delta),
           "raw_len_hist_10nt_both": dict(raw_both),
           "n_span": len(span), "n_span_fail": n_span_fail,
           "identity_median": {"gap50": med(idy_gap), "no_gap": med(idy_nogap)},
           "identity_n": {"gap50": len(idy_gap), "no_gap": len(idy_nogap)},
           "edlib_confirm": {"n": confirm_tot, "confirmed": confirmed},
           "params": {"min_del": a.min_del, "call_del": a.call_del,
                      "primer_fwd": len(FWD), "primer_tail": len(TAIL)}}
    json.dump(out, open(a.out, "w"))
    with open(a.events_out, "w") as fh:
        fh.write("sublib\tgene\tgap_start\tgap_len\tqlen\tmol_len\treflen\n")
        for e in events:
            fh.write("\t".join(map(str, e)) + "\n")
    shutil.rmtree(tmp, ignore_errors=True)
    tot = elig["all"] or 1
    print(f"{a.sublibrary} {os.path.basename(a.fastq)}: {n} subsampled, "
          f"both-primer {ladder.get('both', 0)}, eligible {elig['all']}, "
          f">= {a.call_del} nt deletion {calls[a.call_del]['all']} ({100 * calls[a.call_del]['all'] / tot:.2f}%)")


if __name__ == "__main__":
    main()
