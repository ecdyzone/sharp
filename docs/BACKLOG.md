# Backlog

Ordered by dependency and priority. Items within a tier can be parallelized.

---

## Tier 0 — benchmark priority (team directive: benchmarks first)

- [ ] **Verify MiBIG 4.0 JSON schema**
  Run `python scripts/prepare_mibig_ground_truth.py --inspect data/raw/mibig_json_4.0`.
  Confirm `cluster_id`, `classes`, `taxonomy_name`, and `locus coords` all resolve correctly.
  If any are `None`, edit the `FIELD PATHS` section in `prepare_mibig_ground_truth.py`.
  **Must be done before any real benchmark numbers are trusted.**

- [ ] **`scripts/run_antismash_baseline.py`**
  Run antiSMASH on a genome → parse JSON summary → `data/interim/antismash_predictions.parquet`.
  antiSMASH JSON: `records[].features[]` where `type == "region"`. Set `p_bgc=1.0` (no score).
  Needs antiSMASH installed in the pixi environment (`bioconda::antismash`).

- [ ] **`scripts/run_deepbgc_baseline.py`**
  Run DeepBGC → parse `.bgc.tsv` → `data/interim/deepbgc_predictions.parquet`.
  DeepBGC TSV columns: `sequence_id`, `start`, `end`, `deepbgc_score`, `product_class`.
  Use `deepbgc_score` as `p_bgc`.
  Needs DeepBGC installed (`pip install deepbgc` or check bioconda).

- [ ] **`scripts/prepare_bgcatlas_ground_truth.py`**
  Parse BGC Atlas distribution format → `data/raw/bgcatlas_ground_truth.tsv`.
  Inspect the actual BGC Atlas download format first — it may be GFF, TSV, or JSON.
  Same output schema as `mibig_ground_truth.tsv`. Use as secondary benchmark only.

- [ ] **Literature check: other tools to benchmark against**
  The coworker says "check recent literature if it's worth including more."
  Candidates as of mid-2025: GECCO, SanntiS, BiG-SCAPE (for clustering, not prediction).
  Same pattern: parse output → `predictions.parquet` → `evaluate.py`.

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
  Real tool invocations. Marked `@pytest.mark.slow`. Not run in default `pytest`.
  Start with one test per external tool (Bakta, hmmscan, FIMO, antiSMASH, DeepBGC).

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
- GECCO, SanntiS as additional baselines if literature supports it
- Snakemake / Nextflow DAG for cluster execution

---

## Known tech debt

- `prepare_mibig_ground_truth.py`: schema accessors are written for what we *believe*
  MiBIG 4.0 looks like. Must be verified with `--inspect` on real data before use.
- `extract_embeddings.py`: no resumability. If interrupted on a 50k-protein run,
  restarts from zero.
- `train.py` (not written yet): feature aggregation strategy (mean-pool vs. max-pool
  vs. concat) is an open experiment. Build in a way that's easy to swap.
