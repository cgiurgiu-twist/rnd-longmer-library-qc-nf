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

Agreement is within what a few-hundred-read slice supports: the slice's binomial error on a
6% deletion rate at n≈400 is about ±1.2 points, and every difference is inside that.

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
