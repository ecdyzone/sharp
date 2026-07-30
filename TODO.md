# TODO

## Benchmarks — real data

- [x] Run the full comparison with real data — done for the SCP1 smoke test
      (`AL589148.1`, all three tools). This exposed three defects in the
      benchmark core, now fixed (see `docs/BACKLOG.md` Tier 0). Still to do at
      scale: a genome with more than one coordinate-resolved MiBIG cluster on it.
- [ ] Re-run the comparison on a genome with real recall signal. `AL589148.1`
      carries exactly **one** coordinate-resolved MiBIG cluster, so every recall
      number from it is 0/1 or 1/1. Suggested: *S. coelicolor* A3(2) (AL645882.2)
      — 15 clusters in `streptomyces_ground_truth.tsv`, the most of any contig.
      Always pass `--contigs`, the same file for every tool.
- [ ] S(H)ARP itself can't be benchmarked yet — `predict.py` and the rest of the
      pipeline (`annotate.py` → `train.py`) aren't implemented yet (see CLAUDE.md
      "What is NOT YET IMPLEMENTED").
- [ ] Once real numbers exist, write them up with the MiBIG coordinate-coverage,
      benchmark-scope, and BGC Atlas optimism caveats (CLAUDE.md → "Benchmark
      comparison"). Report `detection` and `reciprocal` recall together — the gap
      between them is the boundary-tightness story (DeepBGC 1.000 vs 0.000 on
      SCP1), and reporting either alone is misleading.
- [ ] Decide whether `min_prediction_frac` should be non-zero for the headline
      table. Currently 0.0, so detection recall ignores how wide a call is; the
      tightness signal lives in `nucleotide.precision` and
      `boundary.median_prediction_coverage` instead. Revisit once a genome with
      more clusters shows whether the two rankings ever disagree.
