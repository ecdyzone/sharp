# %% [markdown]
# # S(H)ARP Benchmarks — Part 1: Ground Truth
#
# **Goal of this notebook:** show how we build the *ground truth* for benchmarking —
# the reference set of "here are real BGCs and where they live" that we score every
# tool against.
#
# The benchmark has **two independent axes**:
#
# 1. **Ground truth** — *what counts as a correct BGC?* → **this notebook**
#    - **MiBIG** — manually curated, our *primary* reference.
#    - **BGC Atlas** — computationally predicted, our *secondary* (noisy) reference.
# 2. **Tools compared** against that ground truth (antiSMASH, DeepBGC, S(H)ARP)
#    → a later notebook.
#
# The key design idea: **both ground-truth sources are normalized to the *same* tiny
# table schema.** Once they look identical, the scoring code (`sharp.evaluate`)
# doesn't care where a cluster came from — MiBIG and BGC Atlas flow through it
# unchanged.
#
# > This is a walkthrough. It shows the real input files and the real output files
# > so you can understand the data without needing the 10 GB dataset on your own
# > machine. For full detail, see the repo's `CLAUDE.md` and `docs/`.

# %% [markdown]
# ## The target schema (what every source becomes)
#
# No matter the source, we produce a **TSV with 5 columns**:
#
# | column | meaning |
# |---|---|
# | `cluster_id` | unique id of the known BGC |
# | `contig` | which sequence it's on |
# | `start`, `end` | position, **0-based half-open** `[start, end)` |
# | `class` | biosynthetic class (PKS, NRPS, terpene, …) |
#
# `sharp.io` gives us a typed `KnownCluster` for exactly this, plus loaders/writers.

# %%
from pathlib import Path

import pandas as pd

from sharp.io import KnownCluster, load_ground_truth_tsv

# All paths are relative to the repo root (run this notebook from notebooks/).
RAW = Path("../data/raw")

# One row of the target schema, as a typed object:
example = KnownCluster(
    cluster_id="BGC0000038",
    contig="AL645882.2",
    start=6890527,
    end=6948414,
    cluster_class="PKS",
)
example

# %% [markdown]
# ---
# ## Source 1 — MiBIG (primary, manually curated)
#
# MiBIG ships as **one JSON file per cluster**. Here's what a raw input entry looks
# like (fields we don't use are trimmed):

# %%
import json

mibig_example = json.loads((RAW / "mibig_json_4.0" / "BGC0000038.json").read_text())

# Show only the fields our parser reads:
{
    "accession": mibig_example["accession"],
    "taxonomy.name": mibig_example["taxonomy"]["name"],
    "biosynthesis.classes": mibig_example["biosynthesis"]["classes"],
    "loci[0]": mibig_example["loci"][0]["accession"],
    "loci[0].location": mibig_example["loci"][0]["location"],
}

# %% [markdown]
# Notice `location: {"from": 6890528, "to": 6948414}`. **MiBIG is 1-based inclusive.**
# Our convention everywhere is **0-based half-open**, so the parser converts on
# ingest (`start - 1`, `end` unchanged). This is the one place that conversion
# happens.
#
# The script that does all this is `scripts/prepare_mibig_ground_truth.py`. Its
# field-reading helpers are importable, so we can show the conversion directly:

# %%
import sys

sys.path.insert(0, "../scripts")
from prepare_mibig_ground_truth import entry_to_clusters  # noqa: E402

# One JSON entry -> one (or more) KnownCluster rows, already coordinate-converted:
entry_to_clusters(mibig_example, "BGC0000038.json")

# %% [markdown]
# → `from: 6890528` (1-based) became `start: 6890527` (0-based). ✅
#
# ### Building the whole MiBIG ground truth
#
# In practice we run the script once over the full dump. Two variants are built:
#
# ```bash
# # All genera
# python scripts/prepare_mibig_ground_truth.py \
#     --input-dir data/raw/mibig_json_4.0 \
#     --output data/raw/mibig_ground_truth.tsv
#
# # Streptomyces only (our organism of interest)
# python scripts/prepare_mibig_ground_truth.py \
#     --input-dir data/raw/mibig_json_4.0 \
#     --output data/raw/streptomyces_ground_truth.tsv \
#     --genus Streptomyces
# ```
#
# Here is the **actual output** on disk:

# %%
mibig_gt = pd.read_csv(RAW / "mibig_ground_truth.tsv", sep="\t")
print(f"MiBIG ground truth: {len(mibig_gt)} clusters (all genera)")
mibig_gt.head()

# %%
strep_gt = pd.read_csv(RAW / "streptomyces_ground_truth.tsv", sep="\t")
print(f"Streptomyces-only ground truth: {len(strep_gt)} clusters")
strep_gt.head()

# %% [markdown]
# ### ⚠️ Important caveat: coordinate coverage
#
# **~53% of *Streptomyces* MiBIG entries have unknown coordinates**
# (`location: {from: 0, to: 0}` — the compound is characterized, but the genomic
# position isn't). A coordinate-based benchmark can't score a cluster with no
# interval, so those are **correctly dropped**.
#
# Consequence: our recall denominator is "recall over *coordinate-resolved* MiBIG,"
# not over all known *Streptomyces* BGCs. This affects **every** tool equally (so the
# *comparison* stays fair), but it caps the absolute numbers. We report it alongside
# every benchmark table.

# %% [markdown]
# ---
# ## Source 2 — BGC Atlas (secondary, computationally predicted)
#
# BGC Atlas ships very differently: **204,661 antiSMASH GenBank (`.gbk`) files**, one
# per predicted region, from metagenome assemblies. Filenames look like:
#
# ```
# MGYA00004361_contig00011.region001.gbk
# └── assembly ──┘└─ contig ─┘└─ region ─┘
# ```
#
# The tricky part: the genomic coordinates are **not** the GenBank `LOCUS` line
# (that's region-local, `1..N`). They live in the antiSMASH structured comment as
# `Orig. start` / `Orig. end` — and, unlike MiBIG, these are **already 0-based
# half-open**, so *no conversion is applied*.
#
# The script `scripts/prepare_bgcatlas_ground_truth.py` isolates all of this. Its
# accessors are importable too — let's read one real file:

# %%
from Bio import SeqIO  # noqa: E402

from prepare_bgcatlas_ground_truth import record_to_cluster  # noqa: E402

gbk_dir = RAW / "complete-bgcs"
example_gbk = sorted(gbk_dir.glob("*.gbk"))[0]
record = SeqIO.read(example_gbk, "genbank")

print("file       :", example_gbk.name)
print("LOCUS coords (region-local, NOT what we want):", f"1..{len(record.seq)}")
print("Orig. coords (genomic, what we DO want)     :",
      record.annotations["structured_comment"]["antiSMASH-Data"]["Orig. start"],
      "->",
      record.annotations["structured_comment"]["antiSMASH-Data"]["Orig. end"])

# One .gbk file -> one KnownCluster row (same type as MiBIG produced):
record_to_cluster(record, example_gbk.name)

# %% [markdown]
# ### Building the BGC Atlas ground truth
#
# Same command shape as MiBIG. The full run streams over all ~204k files; for
# development and demos there's a `--limit` flag so you don't walk 10 GB:
#
# ```bash
# # Full run (~204k rows)
# python scripts/prepare_bgcatlas_ground_truth.py \
#     --input-dir data/raw/complete-bgcs \
#     --output data/raw/bgcatlas_ground_truth.tsv
#
# # Small subset for dev/demo (what produced the file below)
# python scripts/prepare_bgcatlas_ground_truth.py \
#     --input-dir data/raw/complete-bgcs \
#     --output data/raw/bgcatlas_ground_truth.tsv --limit 100
# ```
#
# The **actual output** (100-file sample) on disk:

# %%
bgcatlas_gt = pd.read_csv(RAW / "bgcatlas_ground_truth.tsv", sep="\t")
print(f"BGC Atlas ground truth (sample): {len(bgcatlas_gt)} clusters")
bgcatlas_gt.head()

# %% [markdown]
# ### ⚠️ Important caveat: these labels are predictions
#
# BGC Atlas positives are **antiSMASH predictions, not curated truth**. So "agreeing
# with BGC Atlas" doesn't prove correctness — benchmark numbers against it are
# systematically **optimistic**. We always report it *alongside* MiBIG, never alone.

# %% [markdown]
# ---
# ## Why this matters: one schema → one scorer
#
# Both sources are now the same table. Let's prove it by loading each with the
# **same** loader that the benchmark step uses — `sharp.io.load_ground_truth_tsv` —
# which returns identical `KnownCluster` objects regardless of origin:

# %%
mibig_clusters = load_ground_truth_tsv(RAW / "mibig_ground_truth.tsv")
bgcatlas_clusters = load_ground_truth_tsv(RAW / "bgcatlas_ground_truth.tsv")

print("MiBIG     →", type(mibig_clusters[0]).__name__, ":", mibig_clusters[0])
print("BGC Atlas →", type(bgcatlas_clusters[0]).__name__, ":", bgcatlas_clusters[0])
print()
print("Same type from both sources:",
      type(mibig_clusters[0]) is type(bgcatlas_clusters[0]))

# %% [markdown]
# Because both are lists of `KnownCluster`, the scoring code doesn't branch on the
# source. **Any** tool's predictions get compared to **either** ground truth by the
# exact same function — that's the whole point of normalizing to one schema.
#
# A tiny sanity check with a made-up prediction, using the real metric
# (`reciprocal_overlap`: a prediction "matches" a known cluster if they cover ≥50% of
# each other, on the same contig):

# %%
from sharp.io import PredictedRegion  # noqa: E402
from sharp.metrics import reciprocal_overlap  # noqa: E402

gt = mibig_clusters[0]  # a real known cluster
print("known cluster :", gt)

# A prediction that lands right on top of it → should match:
good = PredictedRegion(
    region_id="pred_1", contig=gt.contig,
    start=gt.start, end=gt.end, p_bgc=0.9,
)
# A prediction far away → should NOT match:
bad = PredictedRegion(
    region_id="pred_2", contig=gt.contig,
    start=gt.start + 10_000_000, end=gt.end + 10_000_000, p_bgc=0.9,
)

print("on-target prediction matches? ", reciprocal_overlap(good, gt, min_frac=0.5))
print("off-target prediction matches?", reciprocal_overlap(bad, gt, min_frac=0.5))

# %% [markdown]
# ---
# ## Recap
#
# - We turn **two very different sources** — MiBIG JSON (curated) and BGC Atlas `.gbk`
#   (predicted) — into **one 5-column table** of known BGCs.
# - The only source-specific work is *reading* each format and *normalizing coordinates*
#   to 0-based half-open (MiBIG needs conversion; BGC Atlas doesn't).
# - Once normalized, `sharp.evaluate` scores **any tool against either ground truth**
#   with the same code.
# - Two caveats travel with the numbers: MiBIG's ~53% coordinate-coverage gap, and
#   BGC Atlas's optimistic (self-predicted) labels.
#
# **Next (Part 2):** run the actual tools (antiSMASH, DeepBGC, S(H)ARP), convert their
# outputs to predictions, and score them against these ground truths.
