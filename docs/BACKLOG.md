# Backlog

Ordered by dependency and priority. Items within a tier can be parallelized.

---

## Tier 0 — benchmark priority (team directive: benchmarks first)

- [x] **Verify MiBIG 4.0 JSON schema** ✅ 2026-07-07
  Ran `--inspect`: `cluster_id`, `classes`, `taxonomy_name`, and locus coords all
  resolve. 1-based→0-based conversion confirmed (e.g. BGC0000002 `from:7362971` →
  `start:7362970`). **Finding:** ~53% of *Streptomyces* entries have `location:
  {from:0,to:0}` (unknown coords) and are correctly dropped → GT is ~430 loci from
  427 clusters, not ~900. Documented as a benchmark caveat in `CLAUDE.md` and
  `PIPELINE.md`. No code change needed; the drop is already logged by the script.

  **Baseline model:** S(H)ARP never invokes antiSMASH/DeepBGC/GECCO (incompatible
  deps → each in its own pixi env via `scripts/setup_<tool>.sh`). The user runs the
  tool; S(H)ARP only parses the output. So each baseline gets a **converter**
  (`convert_<tool>_to_parquet.py`), not a subprocess wrapper. Each isolates column
  names + coordinate base in one place and has an `--inspect` mode; tests use a
  small checked-in real output fixture (no tool execution).

  **Coordinate base — verified 2026-07-15** against a real run of all three tools
  on the same input FASTA (`AL589148.1`), by comparing span (`end-start`) against
  the matching region/cluster `.gbk` LOCUS bp length across every output row:
  antiSMASH **and** DeepBGC are both 0-based half-open (no conversion — refutes
  the old "DeepBGC likely 1-based" guess); GECCO is 1-based inclusive (`start-1`,
  `end` unchanged — confirms the old guess). Full evidence in `CLAUDE.md` →
  "Baseline integration". Converter plans below are finalized against real output,
  not yet implemented.

- [x] **`scripts/convert_antismash_to_parquet.py`** ✅ 2026-07-15
  Parses `sequence.json` → `data/interim/antismash_predictions.parquet`.
  `records[].features[]` where `type == "region"`. `location` is a string
  `"[start:end](strand)"` (regex parse, not plain ints), 0-based half-open, no
  conversion. `region_id = f"{contig}.region{region_number:03d}"` (matches the
  `.gbk` filename; `region_number` alone resets per contig). `predicted_class` =
  `;`-joined `qualifiers['product']` (hybrid regions list multiple products).
  `p_bgc=1.0` (no score). Loops all `records`, not just the first (multi-contig).
  `--inspect` mode; tests: `tests/test_convert_antismash.py` (23 tests) against a
  checked-in trimmed real fixture (`tests/fixtures/antismash_sequence.json`) and
  verified live against the full real run output.

- [x] **`scripts/convert_deepbgc_to_parquet.py`** ✅ 2026-07-15
  Parses `<prefix>.bgc.tsv` → `data/interim/deepbgc_predictions.parquet`.
  Coordinate columns are `nucl_start`/`nucl_end` (not `start`/`end`), 0-based
  half-open, no conversion. `region_id = bgc_candidate_id` (already unique).
  `p_bgc = deepbgc_score`. `predicted_class = product_class` — confirmed to
  exist (28-col real header) but frequently empty string (4/5 rows in the
  verification run) → mapped to `None`. `--inspect` mode; tests:
  `tests/test_convert_deepbgc.py` (25 tests) against a checked-in real fixture
  (`tests/fixtures/deepbgc_out.bgc.tsv`, unmodified, 5 rows) and verified live
  against the full real output.

- [ ] **`scripts/convert_gecco_to_parquet.py`** — plan finalized 2026-07-15
  Parse `<genome>.clusters.tsv` → `data/interim/gecco_predictions.parquet`.
  Columns confirmed: `sequence_id`, `cluster_id`, `start`, `end`, `average_p`,
  `max_p`, `type`, plus per-class probability columns. Coordinates are 1-based
  inclusive — convert `start-1`, `end` unchanged. `region_id = cluster_id`
  (already unique). `p_bgc = average_p`. `predicted_class = type`, kept as-is even
  though it was `"Unknown"` for 5/5 rows in the verification run (argmax over the
  per-class probability columns is a possible v2 improvement, not done now).

- [x] **`scripts/prepare_bgcatlas_ground_truth.py`** ✅ 2026-07-07
  Parses the BGC Atlas `complete-bgcs` dump (204,661 antiSMASH `.gbk` files, one
  region per file) → `data/raw/bgcatlas_ground_truth.tsv` (same schema as MiBIG GT).
  Verified schema on real data: genomic coords come from the antiSMASH `Orig.
  start`/`Orig. end` structured-comment fields (NOT the region-local LOCUS coords),
  and are already 0-based half-open (`end-start == len(seq)`), so no conversion.
  `cluster_id` = filename stem (unique, includes `.regionNNN`); `contig` =
  `<MGYA assembly>_<rec.id>` (assembly-qualified, since rec.id repeats across
  assemblies). `--limit N` for dev/tests; `--inspect` to re-verify schema. Tests:
  `tests/test_prepare_bgcatlas.py` (24, synthetic gbk fixtures). Secondary/noisy GT.

- [ ] **Literature check: other tools to benchmark against**
  The coworker says "check recent literature if it's worth including more."
  antiSMASH, DeepBGC, GECCO are already in scope (converters above). Remaining
  candidates as of mid-2025: SanntiS, BiG-SCAPE (for clustering, not prediction).
  Same pattern: run tool in its own env → convert output → `evaluate.py`.

- [ ] **`scripts/build_sarp_hmm.py`**
  Collect SARP seed sequences → align → `hmmbuild` → `data/raw/sarp_models.hmm`.
  Blocks step 02a of the main pipeline (not the benchmark).

---

## Tier 1 — core pipeline steps (implement in order)

- [ ] **`sharp/annotate.py`** — Bakta wrapper
  Inputs: `genome.fasta`
  Outputs: `proteins.faa`, `annotated.gbk`, `genes.gff`
  See `CLAUDE.md` step 1 for full spec.

- [ ] **`sharp/detect_sarp.py`** — hmmscan wrapper + tblout parser
  Inputs: `proteins.faa`, `sarp_models.hmm`
  Output: `anchors_sarp.tsv`

- [ ] **`sharp/detect_heptarepeats.py`** — FIMO wrapper + output parser
  Inputs: `genome.fasta`, motif file
  Output: `anchors_heptarepeats.tsv`
  **Caution:** verify whether FIMO output coords are 0-based or 1-based before
  writing the parser. Use `--inspect`-style debug output.

- [ ] **`sharp/merge_anchors.py`** — pure concat + dedup
  Inputs: `anchors_sarp.tsv`, `anchors_heptarepeats.tsv`
  Output: `anchors.tsv`

- [ ] **`sharp/extract_neighborhood.py`** — BioPython GenBank parsing
  Inputs: `anchors.tsv`, `annotated.gbk`
  Outputs: `neighborhoods.tsv`, `neighborhood_proteins.faa`, `neighborhood_dna.fna`
  Key correctness concern: the FASTA headers must be `>PROTEIN_ID region_id=R001`
  exactly, or `parse_fasta` in `io.py` will skip them.

- [ ] **`sharp/annotate_domains.py`** — hmmscan domtblout parser
  Inputs: `neighborhood_proteins.faa`, `pfam_models.hmm`
  Output: `domains.tsv`

---

## Tier 2 — ML steps (can start in parallel with Tier 1 using mock data)

- [ ] **`sharp/extract_kg_features.py`** ⭐
  Inputs: `neighborhoods.tsv`, `domains.tsv`, `kg.gpickle`
  Output: `kg_features.parquet`
  Can be developed with mock `neighborhoods.tsv` and `domains.tsv` before
  Tier 1 is done. Needs `kg.gpickle` from `build_kg.py`.

- [ ] **`scripts/build_kg.py`** ⭐
  Inputs: `mibig_json_4.0/`, BGC Atlas data
  Output: `kg.gpickle`
  Unblocks `extract_kg_features.py`.

- [ ] **`sharp/train.py`** ⭐
  Inputs: `embeddings.parquet`, `domains.tsv`, `kg_features.parquet`, `mibig_ground_truth.tsv`
  Outputs: `model.pkl`, `metrics.json`, `feature_importance.tsv`
  Can prototype feature aggregation logic with mock data before Tier 1 is done.

- [ ] **`sharp/predict.py`**
  Inputs: `model.pkl` + same features as train
  Output: `predictions.parquet`
  Depends on `train.py` being done.

---

## Tier 3 — output + post-processing

- [ ] **`sharp/filter.py`**
  Heuristic post-filter on predictions.

- [ ] **`sharp/generate_report.py`**
  Jinja2 HTML report.

---

## Tier 4 — validation infrastructure

- [ ] **`tests/integration/`** directory with slow tests
  Real invocations of the **pipeline's own** tools (Bakta, hmmscan, FIMO), which
  live in the S(H)ARP env. Marked `@pytest.mark.slow`. Not run in default `pytest`.
  NB: the baseline tools (antiSMASH/DeepBGC/GECCO) are *not* tested here — they run
  in separate envs; their converters are tested against checked-in output fixtures.

- [ ] **End-to-end smoke test on a known Streptomyces genome**
  Suggested: *S. coelicolor* A3(2) (NCBI: AL645882). MiBIG has ~20 clusters from it.
  Run full pipeline + all baselines → compare three `benchmark.json` files.
  This is the first real benchmark table for the paper.

---

## Future / not yet scoped

- BGC Atlas as training data augmentation (noisy positives for the classifier)
- Per-class benchmark breakdown — "which BGC class is each tool best at?"
- Evo (nucleotide language model) embeddings — `neighborhood_dna.fna` already generated
- GNN embeddings from knowledge graph
- Multi-class classification
- SanntiS as an additional baseline if literature supports it
- Snakemake / Nextflow DAG for cluster execution

---

## Known tech debt

- `prepare_mibig_ground_truth.py`: schema accessors are written for what we *believe*
  MiBIG 4.0 looks like. Must be verified with `--inspect` on real data before use.
- `extract_embeddings.py`: no resumability. If interrupted on a 50k-protein run,
  restarts from zero.
- `train.py` (not written yet): feature aggregation strategy (mean-pool vs. max-pool
  vs. concat) is an open experiment. Build in a way that's easy to swap.
