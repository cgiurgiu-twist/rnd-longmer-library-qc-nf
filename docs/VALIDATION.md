# Validation against the Databricks analysis

The reference implementation of this QC is the Databricks Spark analysis used for the
AILK013 and AILK014 750mer deliveries. This pipeline is a port of that analysis, and the
port was checked against it on the same sequencing run.

## Run used

Q-724749, AILK014-007 / -008 / -009, three PromethION flowcells (FLO-PRO114M, SQK-LSK114,
ONT barcoding off), 93,469,486 reads, 68,189 designed constructs.

## Comparison

The Databricks numbers are full depth against the whole reference. The pipeline numbers are
from `-profile local`: a 720-construct slice of the same reference (six high-density wells,
two per sub-library) and the ~1,800 reads out of 152,000 scanned that map to it.

| Metric | Databricks, full depth | Pipeline, test slice |
|---|---|---|
| Correctly paired, -007 | 99.46% | 99.86% |
| Correctly paired, -008 | 99.47% | 100.00% |
| Correctly paired, -009 | 99.70% | 99.72% |
| ≥50 nt block deletion, -007 | 6.14% | 6.79% |
| ≥50 nt block deletion, -008 | 6.98% | 7.48% |
| ≥50 nt block deletion, -009 | 5.88% | 5.34% |
| Intact full-length (both-primer), -007 | 93.8% | 93.2% |
| Intact full-length (both-primer), -008 | 93.0% | 92.5% |
| Intact full-length (both-primer), -009 | 94.1% | 94.7% |

Largest differences: pairing 0.53 pp (-008), deletion rate 0.65 pp (-007). Both are inside
what a slice this small supports — on a 6% rate at n≈400 the binomial standard error is
about 1.2 pp, so the 95% interval is roughly ±2.3 pp.

The slice is 1,812 reads in total and -009 contributes only 364, so its "99.72% correctly
paired" is a single mispaired read. Read this table as "the ported callers behave the same
on the same reads", not as a precision measurement.

## What this does and does not establish

It establishes that the ported classifier and deletion caller reproduce the reference
implementation's rates on the same reads — the half-split rule, the margin test, the
class definitions, the both-primer filter and the splice-gap call all behave the same.

It does **not** establish agreement on uniformity, dropout, 95/5 or screening depth. Those
are depth-dependent and a slice cannot speak to them; on the slice most variants have a
handful of reads and the 95/5 ratio is undefined. Those metrics were verified by
inspection of the merge arithmetic against the Databricks aggregation, not by a numeric
comparison.

It is also not a test of the classifier's correctness in absolute terms — that is what the
per-run controls in `metrics/controls.json` are for, and they run on every execution.

## Re-running the comparison

```bash
nextflow run . -profile local
python3 - <<'PY'
import json
d = json.load(open('results_test/metrics/metrics_core.json'))
s = json.load(open('results_test/metrics/blockdel_scope.json'))
print({o: round(v['clean_designed_pct_of_assigned'], 2) for o, v in d['per_sublib'].items()})
print({o: round(v['pct'], 2) for o, v in s['per_sublib_any_cls'].items()})
print({o: round(v, 1) for o, v in s['fl_readlen_nogap_by_sublib_both'].items()})
PY
```

Expected, for the bundled fixtures:

```
{'AILK014-007': 99.86, 'AILK014-008': 100.0, 'AILK014-009': 99.72}
{'AILK014-007': 6.79, 'AILK014-008': 7.48, 'AILK014-009': 5.34}
{'AILK014-007': 93.2, 'AILK014-008': 92.5, 'AILK014-009': 94.7}
```

A change to the classifier, the half-split rule or the deletion caller that moves these is
a real behaviour change and needs to be explained in the PR, not re-baselined silently.

## Reference-builder verification (design side)

`bin/make_reference.py` builds the reference CSV from the delivered order files and checks it
against the MOP manifests and the writer files that went to the printer. On the Q-724749
build (AILK014-007/-008/-009) it reports:

| sub-library | manifest well agreement | writer rows identical to design |
|---|---|---|
| AILK014-007 | 22,391 / 22,391 | 44,782 / 44,782 |
| AILK014-008 | 23,000 / 23,000 | 46,000 / 46,000 |
| AILK014-009 | 22,798 / 22,798 | 45,596 / 45,596 |

Writer files carry two replicate rows per oligo, hence ~2× the construct count, and encode
every `A` as `8`; that substitution is reversed before comparison.

The manifest records each check **per sub-library** with `run`, `n_checked`, `n_ok` /
`n_identical`, `n_mismatched` and `n_unrecognised_names` — not a flat list of check names. A
check whose input file is absent is recorded as `"run": false` with a reason, and a check that
reads its file but matches **zero** designed oligos is recorded as a problem rather than a
pass. That guard exists because an earlier version keyed the design lookup on `Gene` while the
writer files key on `Oligo_ID`: it matched nothing, reported `n_checked: 0, n_mismatched: 0`,
and read as a clean result. Zero mismatches out of zero comparisons is not evidence.
