#!/usr/bin/env python3
"""merge_events.py — concatenate the per-shard deletion-event samples into one table."""
import argparse, csv, glob, random, sys

ap = argparse.ArgumentParser()
ap.add_argument("--events", nargs="+", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--max", type=int, default=300000, help="reservoir-sample down to this many events")
ap.add_argument("--seed", type=int, default=20260927)
a = ap.parse_args()

rng = random.Random(a.seed)
keep, n = [], 0
for f in sorted(set(a.events)):
    with open(f) as fh:
        next(fh, None)
        for line in fh:
            n += 1
            if len(keep) < a.max:
                keep.append(line)
            else:
                j = rng.randrange(n)
                if j < a.max:
                    keep[j] = line
with open(a.out, "w") as fh:
    fh.write("sublib,gene,gap_start,gap_len,qlen,reflen\n")
    for l in keep:
        fh.write(",".join(l.rstrip("\n").split("\t")) + "\n")
print(f"{n} events seen, {len(keep)} written to {a.out}")
