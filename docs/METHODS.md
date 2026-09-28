# Methods — what each number means

## Read classification

Each read is aligned once (minimap2 `map-ont`, via `mappy`) against a combined index that
holds, for every designed construct, a **VH-side half** and a **VL-side half**. The halves
are cut at the midpoint of the shared intermediate:

```
5' flank (33) | VH | intermediate (162) | VL | 3' flank (33)
              └────── H half ──────┘└────── L half ──────┘
                        cut at 33 + vh_len_nt + 81
```

The cut sits on the crossover hotspot: in a printed molecule a VH(A)+VL(B) read can only
arise from a template switch somewhere between the variable regions, and the intermediate
is the stretch in between. It is constant at the protein level but codon-diversified per
gene, so it still discriminates.

A half is assigned only if its best hit beats the runner-up by `--margin` (15) matching
bases over an alignment of at least `--min_mlen` (120) nt.

| class | meaning |
|---|---|
| `clean_designed` | both halves are the same designed gene |
| `within_well` | different genes, same high-density well → recombination inside the well |
| `cross_well` | different wells, same sub-library |
| `cross_sublib` | different sub-libraries — with one flowcell each, this is a control, not biology |
| `one_half` | only one half assignable |
| `unmapped` | neither half assignable |

**The margin is load-bearing and its cost is not neutral.** A half carrying an internal
recombination is exactly the kind of half that matches two references closely and fails the
margin test, so the excluded reads are plausibly enriched for recombinants and the reported
chimera rate is likely conservative.

## Why not soft-clips

In a conserved-framework antibody library a chimera aligns end-to-end to one designed gene
with no soft clip and no supplementary alignment, so an aligner-rescue caller scores it
clean. Measured undercounts on other Twist builds range from ~9× to ~80×. Content
deconvolution is not an optimisation here — it is the only method that sees the defect.

## Internal block deletions

A read must prove it still has both ends before anything is said about its middle:

1. `cutadapt` with a linked adapter requiring **both** primers, `--action=retain` so the
   read still spans the design. Reads that fail are end-truncated and are profiled
   separately (Section 8 of the report), never counted as deletions.
2. `minimap2` splice preset against the read's **own** assigned construct. A pilot that
   took the best hit across all constructs called a gap in 13% of reads at 0.75 identity —
   those were reads landing on a near-neighbour variant with the splice model bridging CDR
   mismatches as false introns. Per-construct alignment removes that failure mode.
3. A call is one internal gap ≥ `--call_del` (50 nt), cross-checked by an `edlib` global
   alignment with no splice model.

50 nt is where the controls give full sensitivity and zero false calls. `edlib` on its own
is not a caller — it produces false gaps at ≥3% read error — it is only a confirmation.

## Full length

Four definitions, three denominators, because they disagree by more than 10 points on the
same reads:

| definition | what it misses |
|---|---|
| alignment span ≥90% of design | **straddles an internal deletion — cannot see this defect at all** |
| read length ≥90% of design | tolerates the loss of up to ~68 nt on a 687 nt construct |
| no internal deletion ≥50 nt | says nothing about the ends |
| read length ≥90% **and** no deletion ≥50 nt | the one we quote |

Quoted against reads that retain both primer sites — the only population in which an
internal loss is scorable.

## Uniformity, dropout and screening depth

Per-variant counts come from `clean_designed` reads only. CV includes dropped variants as
zeros. The 95/5 ratio is **undefined** when the 5th percentile is zero, and is reported as
`n/a` with that reason rather than as a large number.

Recovery is Poisson: `coverage(N) = mean_i [1 − exp(−N · p_i · f_i)]`, computed per variant
from its own abundance and its own measured deletion rate. Because the intact fraction is
an upper bound, the cell counts are **minima**. A target that cannot be reached — too many
variants with no reads — returns `n/a`, not the search ceiling.

## Calibration

Run on every execution, against this library's own constructs:

- **Classifier negatives**: simulated clean reads at 1/2/3/5/8% error, uniform and
  clustered. Clustered matters: nanopore error arrives in bursts and a caller tuned only
  against uniform error will manufacture events in a burst. Zero calls gives a 95% upper
  bound of 3/n.
- **Classifier positives**: injected within-well VH/VL swaps.
- **Deletion caller**: the both-primer filter against truncated and intact reads; error-only
  reads for false calls; programmed deletions of known size and position for sensitivity
  and placement.

A library with more well-mates, or shorter constructs, is a harder problem — that is why
the operating point is re-measured per run rather than inherited from a previous one.

## Known limitations

- **Intra-segment mosaics are not measured.** The Databricks analysis includes a window-vote
  breakpoint detector within a half; it is not ported here. Intact fractions from this
  pipeline are therefore upper bounds that do not account for within-half recombination.
- **`cross_well` is reported, not adjudicated.** On these builds most cross-well reads are
  ligation-fused molecules from sequencing prep rather than chimeric constructs; the
  read-length-by-class table is what distinguishes them, and the call is left to the reader.
- **Rarefaction is modelled, not re-counted** — Poisson thinning of the observed per-variant
  counts, which is exact for uniform read subsampling but is not a second pass over the data.
