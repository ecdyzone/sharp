# Pipeline — S(H)ARP

Full data flow from genome to BGC predictions. Each numbered step maps to
a module in `src/sharp/`. Status: ✅ done · 🔲 not started.

---

## Pre-requisites (one-time setup, not part of the main run)

These artifacts are built once and versioned. They're inputs to the pipeline,
not produced by it.

| Artifact | How built | Used by |
|---|---|---|
| `data/raw/sarp_models.hmm` | `scripts/build_sarp_hmm.py` (hmmbuild on SARP seed sequences) | `detect_sarp.py` |
| `data/raw/pfam_models.hmm` | Download from InterPro/Pfam | `annotate_domains.py` |
| `data/raw/kg.gpickle` | `scripts/build_kg.py` (NetworkX, from MiBIG + BGC Atlas) | `extract_kg_features.py` |
| `data/raw/mibig_ground_truth.tsv` | `scripts/prepare_mibig_ground_truth.py` | `train.py`, `evaluate.py` |
| `data/raw/mibig_json_4.0/` | `scripts/download_mibig.sh` | `prepare_mibig_ground_truth.py` |
| `data/raw/mibig_prot_seqs_4.0.fasta` | `scripts/download_mibig.sh` | `build_kg.py`, future |

The **knowledge graph** (`kg.gpickle`) encodes:
- Nodes: BGC clusters (from MiBIG), domain types (from Pfam), biosynthetic classes
- Edges: "cluster has domain", "cluster has class", "domain co-occurs with domain"
- Built with NetworkX, serialized with `pickle` (`gpickle` = gzipped pickle)
- `scripts/build_kg.py` is not yet written — see backlog

---

## Step 00 — Input

```
data/raw/genome.fasta
```

Raw genome of an actinomycete in FASTA format. No annotation assumed.
May be multi-contig. Obtained from NCBI or own sequencing.

---

## Step 01 — Genome annotation ✅ (designed, 🔲 implement)

```
genome.fasta
    └─► annotate.py (Bakta)
            ├─► proteins.faa         data/interim/
            ├─► annotated.gbk        data/interim/
            └─► genes.gff            data/interim/
```

Bakta detects ORFs, translates proteins, assigns functional annotations.
`annotated.gbk` is the canonical inter-step format (coordinates + sequences).
`genes.gff` is auxiliary (IGV, external tools) — not consumed downstream.

---

## Step 02 — Anchor detection (two parallel frentes) 🔲

### 02a — SARP by protein HMM

```
proteins.faa + sarp_models.hmm
    └─► detect_sarp.py (hmmscan)
            └─► anchors_sarp.tsv     data/interim/
```

Columns: `protein_id, contig, start, end, strand, score, type=SARP`
One row per significant HMM hit. E-value threshold: 1e-5 (configurable).

Three SARP architectures targeted:
- Small: HTH-BTAD only
- Medium: HTH-BTAD + NB-ARC
- Large: HTH-BTAD + NB-ARC + TPR / AAA / LuxR

### 02b — Heptarepeats by DNA motif

```
genome.fasta + afsR-box PWM
    └─► detect_heptarepeats.py (FIMO)
            └─► anchors_heptarepeats.tsv    data/interim/
```

Columns: same schema as 02a, `type=heptarepeat`
FIMO p-value threshold: 1e-4 (configurable).

### 02c — Merge

```
anchors_sarp.tsv + anchors_heptarepeats.tsv
    └─► merge_anchors.py
            └─► anchors.tsv          data/interim/
```

Concatenate, deduplicate, sort. Pure function.

---

## Step 03 — Neighborhood extraction 🔲

```
anchors.tsv + annotated.gbk
    └─► extract_neighborhood.py
            ├─► neighborhoods.tsv          data/interim/
            ├─► neighborhood_proteins.faa  data/interim/
            └─► neighborhood_dna.fna       data/interim/
```

Window: ±20 genes **or** ±20 kb (whichever is larger).
Merge anchors within 50 kb → one region (avoids double-counting).

`neighborhoods.tsv` columns: `region_id, contig, start, end, anchor_ids, n_proteins`
FASTA headers: `>PROTEIN_ID region_id=R001` (required by `parse_fasta` in io.py)
`neighborhood_dna.fna` — reserved for Evo (not consumed by MVP).

---

## Step 04 — Domain annotation 🔲

```
neighborhood_proteins.faa + pfam_models.hmm
    └─► annotate_domains.py (hmmscan)
            └─► domains.tsv          data/interim/
```

Columns: `protein_id, region_id, domain, e_value, start, end`
One row per domain hit per protein. `region_id` recovered from FASTA header.
Pfam domtblout format parsed.

---

## Step 05 — Embedding extraction ✅

```
neighborhood_proteins.faa
    └─► extract_embeddings.py (ESM-2)
            └─► embeddings.parquet   data/interim/
```

Columns: `protein_id, region_id, embedding[320]` (320-dim for default 8M model; 1280-dim for 650M)
Mean-pool over residue tokens. CLS/EOS/pad masked out of the mean.
Model: `esm2_t6_8M_UR50D` (default, laptop-friendly).

---

## Step 06 — Knowledge graph features ✅ (designed, 🔲 implement)

```
neighborhoods.tsv + domains.tsv + kg.gpickle
    └─► extract_kg_features.py
            └─► kg_features.parquet  data/interim/
```

Columns: `region_id, n_similar_clusters, modal_class, has_large_sarp, ...`
Already aggregated per region (unlike embeddings and domains which are per protein).

For each region: query KG for clusters sharing domain architecture → count matches,
extract modal class, flag presence of large SARP (HTH-BTAD + ≥2 extra domains).

---

## Step 07 — Training ✅ (designed, 🔲 implement)

```
embeddings.parquet + domains.tsv + kg_features.parquet + mibig_ground_truth.tsv
    └─► train.py (LightGBM)
            ├─► model.pkl                data/processed/
            ├─► metrics.json             data/processed/
            └─► feature_importance.tsv   data/processed/
```

Feature aggregation (proteins → regions):
1. `embeddings.parquet` → mean-pool by `region_id` → vector per region
2. `domains.tsv` → one-hot by `region_id` → sparse vector per region
3. `kg_features.parquet` → already per region, join directly

Labels: region is positive if it overlaps a MiBIG cluster with ≥50% reciprocal
overlap (uses `reciprocal_overlap` from `metrics.py`).

Binary classifier (MVP). k-fold CV (k=5). Output: probability + metrics.

---

## Step 08 — Prediction 🔲

```
model.pkl + (same features as train)
    └─► predict.py
            └─► predictions.parquet  data/interim/
```

Columns: `region_id, contig, start, end, p_bgc, predicted_class`
Must match `PredictedRegion` schema in `io.py`.

---

## Step 09 — Filtering 🔲

```
predictions.parquet + domains.tsv
    └─► filter.py
            └─► filtered_predictions.parquet   data/interim/
```

Heuristic rules (MVP):
- Drop if `p_bgc < 0.5`
- Drop if region contains only ribosomal domains
- Drop if `n_proteins < 3`

---

## Step 10 — Report 🔲

```
filtered_predictions.parquet + neighborhoods.tsv + domains.tsv
    └─► generate_report.py (Jinja2)
            └─► report.html          data/processed/
```

One section per BGC: coordinates, class, score, domain architecture.

---

## Benchmark loop ✅ (parallel to main pipeline, team priority)

The same `evaluate.py` handles all three tools — each just needs its output
converted to `predictions.parquet` first.

```
                         mibig_ground_truth.tsv   (primary)
                         bgcatlas_ground_truth.tsv (secondary, noisy)
                                │
          ┌─────────────────────┼─────────────────────┐
          │                     │                     │
antismash_predictions.parquet   │         deepbgc_predictions.parquet
          │             predictions.parquet           │
          │               (S(H)ARP)                   │
          └──────────────────── ┼ ────────────────────┘
                                │
                           evaluate.py
                                │
               ┌────────────────┼────────────────┐
               │                │                │
    benchmark_antismash.json  benchmark_sharp.json  benchmark_deepbgc.json
```

Conversion scripts (not yet written, Tier 0 priority):
- `scripts/run_antismash_baseline.py` — antiSMASH JSON → predictions.parquet
- `scripts/run_deepbgc_baseline.py` — DeepBGC TSV → predictions.parquet
- `scripts/prepare_bgcatlas_ground_truth.py` — BGC Atlas → ground_truth.tsv

BGC Atlas note: computationally predicted, no manual curation. Benchmark numbers
against it are systematically optimistic. Always report alongside MiBIG numbers.

See `docs/ARCHITECTURE.md#metrics` for the methodological choices (match definition,
recovery semantics, threshold recording).
