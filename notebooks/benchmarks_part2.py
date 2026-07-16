# %% [markdown]
# # S(H)ARP Benchmarks — Part 2: Tool Outputs → predictions.parquet
#
# **Goal of this notebook:** show how raw output from the three baseline tools
# (antiSMASH, DeepBGC, GECCO) gets normalized into predictions we can score.
#
# Recall the benchmark's two independent axes (from Part 1):
#
# 1. **Ground truth** — *what counts as a correct BGC?* (MiBIG / BGC Atlas) → Part 1.
# 2. **Tools compared** against that ground truth → **this notebook.**
#
# Same design idea as Part 1, mirrored: **every tool's output is normalized to the
# same tiny table schema** (`PredictedRegion`). Once they look identical, `sharp.evaluate`
# doesn't care which tool produced a prediction — antiSMASH, DeepBGC, GECCO, and
# S(H)ARP itself all flow through it unchanged.
#
# S(H)ARP **never runs these tools**. Each has mutually incompatible dependencies
# and installs into its own isolated pixi env (`scripts/setup_<tool>.sh`). You run
# the tool yourself — its own env, HPC, or a container — and S(H)ARP only *parses
# the output files it leaves behind*. Each tool gets one **converter** script, not
# a subprocess wrapper.
#
# All three tools were run once, for real, on the same input: `AL589148.1`, the
# SCP1 plasmid of *Streptomyces coelicolor* — small, well-characterized, a good
# smoke-test genome. The files loaded below are real (sometimes trimmed) tool
# output checked into `tests/fixtures/`, not synthetic data.

# %% [markdown]
# ## The target schema (what every tool's output becomes)
#
# | column | meaning |
# |---|---|
# | `region_id` | unique id of the predicted region |
# | `contig` | which sequence it's on |
# | `start`, `end` | position, **0-based half-open** `[start, end)` |
# | `p_bgc` | the tool's BGC probability/confidence |
# | `predicted_class` | biosynthetic class, if the tool provides one |
#
# `sharp.io` gives us a typed `PredictedRegion` for exactly this, plus loaders/writers
# — the same role `KnownCluster` played for ground truth in Part 1.

# %%
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, "../scripts")
from sharp.io import PredictedRegion, load_predictions_parquet, write_predictions_parquet

# Fixtures are real (or trimmed-real) tool output checked into the repo, so this
# notebook runs anywhere without needing a 10 GB external benchmark directory.
FIXTURES = Path("../tests/fixtures")

# One row of the target schema, as a typed object:
example = PredictedRegion(
    region_id="AL589148.1.region001",
    contig="AL589148.1",
    start=201195,
    end=222794,
    p_bgc=1.0,
    predicted_class="terpene",
)
example

# %% [markdown]
# ---
# ## Detour: verifying coordinate conventions before trusting any parser
#
# Before writing a single converter, we had to answer one question per tool:
# *is `start`/`end` 0-based half-open (our convention) or 1-based inclusive
# (MiBIG's convention, see Part 1)?* Getting this wrong silently shifts every
# prediction by one base pair — the kind of bug that doesn't crash, it just
# quietly corrupts every overlap calculation.
#
# We didn't assume — we measured. Each tool also writes a GenBank (`.gbk`) file
# per region/cluster, and a `.gbk`'s `LOCUS` line reports its exact length. So for
# every row, we can compare:
#
# - `span = end - start` (from the tool's own table/JSON)
# - `LOCUS bp` (from the matching `.gbk` file)
#
# If `span == LOCUS bp`, the source is **0-based half-open** (our convention,
# no conversion needed). If `span == LOCUS bp - 1`, the source is **1-based
# inclusive** (needs `start - 1` on ingest, same as MiBIG in Part 1).
#
# > **Caution that bit us mid-investigation:** checking only *one* row can look
# > consistent with either convention if you eyeball the wrong direction of the
# > off-by-one. We caught this on a re-check across GECCO's full 5-row output —
# > always verify the span-vs-LOCUS relationship across *every* row, not just one.

# %%
# GECCO cluster_1: real numbers from sequence.clusters.tsv vs. the matching .gbk
gecco_start, gecco_end = 20274, 53842
gecco_locus_bp = 33569

span = gecco_end - gecco_start
print(f"span (end - start)      = {span}")
print(f"LOCUS bp                = {gecco_locus_bp}")
print(f"span == LOCUS bp?       {span == gecco_locus_bp}   (would mean 0-based half-open)")
print(f"span == LOCUS bp - 1?   {span == gecco_locus_bp - 1}   (means 1-based inclusive)")

# %% [markdown]
# `33568 == 33569 - 1` → GECCO's `start`/`end` are **1-based inclusive**. This
# check, repeated across all three tools and every row (not just one), gave us:
#
# | Tool | Coordinate base | Conversion |
# |---|---|---|
# | antiSMASH | 0-based half-open (verified) | none |
# | DeepBGC | 0-based half-open (verified) | none — *refutes* our initial guess |
# | GECCO | 1-based inclusive (verified) | `start - 1`, `end` unchanged — *confirms* our initial guess |
#
# One hypothesis held, one didn't — which is exactly why we measured instead of
# assuming. Full evidence lives in `CLAUDE.md` → "Baseline integration".

# %% [markdown]
# ---
# ## Source 1 — antiSMASH
#
# antiSMASH writes a JSON summary (named after the input FASTA, e.g.
# `sequence.json`), plus one GenBank file per predicted region. The summary's
# `records[].features[]` list mixes every feature type (genes, domains, ...) —
# we filter for `type == "region"`:

# %%
import json

antismash_data = json.loads((FIXTURES / "antismash_sequence.json").read_text())
region_features = [
    f for f in antismash_data["records"][0]["features"] if f["type"] == "region"
]
region_features[0]

# %% [markdown]
# Two things to notice:
#
# - `location` is a **string** — `"[201195:222794](+)"` — not plain integers.
#   Needs a regex, not `int(feature["location"])`.
# - `qualifiers.region_number` (`["1"]`) resets **per contig**, so it can't be
#   used as a global id on its own. We build `region_id` as
#   `f"{contig}.region{region_number:03d}"`, which conveniently matches the
#   tool's own `.gbk` filename convention (`AL589148.1.region001.gbk`) — so a
#   row can always be traced back to its source file.
#
# The converter's own parsing function does exactly this:

# %%
sys.path.insert(0, "../scripts")
from convert_antismash_to_parquet import record_to_regions  # noqa: E402

antismash_regions = []
for record in antismash_data["records"]:
    antismash_regions.extend(record_to_regions(record))
antismash_regions

# %% [markdown]
# The second region here (`region002`) is a hybrid — antiSMASH listed two
# products (`furan`, `butyrolactone`). We join them with `;` rather than
# picking one and silently dropping information.
#
# antiSMASH gives no probability score (it's rule-based), so `p_bgc` is fixed
# at `1.0` for every row.
#
# ### Running it for real
#
# ```bash
# python scripts/convert_antismash_to_parquet.py --inspect <antismash output dir>
# python scripts/convert_antismash_to_parquet.py \
#     --input <antismash output dir> --output data/interim/antismash_predictions.parquet
# ```

# %% [markdown]
# ---
# ## Source 2 — DeepBGC
#
# DeepBGC writes a plain TSV, `<prefix>.bgc.tsv` — one row per candidate BGC,
# no GenBank parsing needed. The real header has 28 columns; only a few matter:

# %%
deepbgc_df = pd.read_csv(FIXTURES / "deepbgc_out.bgc.tsv", sep="\t")
print(f"{len(deepbgc_df)} candidate rows, {len(deepbgc_df.columns)} columns")
deepbgc_df[["sequence_id", "bgc_candidate_id", "nucl_start", "nucl_end",
            "deepbgc_score", "product_class"]]

# %% [markdown]
# Two things to notice:
#
# - The coordinate columns are `nucl_start`/`nucl_end` — **not** `start`/`end`
#   (unlike GECCO below, which uses those exact names for a *different*
#   convention). Column names alone don't tell you the convention; that's why
#   we measured against `.gbk` LOCUS lengths instead of assuming.
# - `product_class` is **frequently blank** — 4 of these 5 real rows have no
#   confident class call. We map blank to `None`, not to an error; downstream
#   code must expect a missing class as the normal case, not the exception.

# %%
from convert_deepbgc_to_parquet import rows_to_regions as deepbgc_rows_to_regions  # noqa: E402

with open(FIXTURES / "deepbgc_out.bgc.tsv", newline="") as fh:
    import csv
    deepbgc_rows = list(csv.DictReader(fh, delimiter="\t"))
deepbgc_regions = deepbgc_rows_to_regions(deepbgc_rows)
deepbgc_regions

# %% [markdown]
# ### Running it for real
#
# ```bash
# python scripts/convert_deepbgc_to_parquet.py --inspect <deepbgc output dir>
# python scripts/convert_deepbgc_to_parquet.py \
#     --input <deepbgc output dir> --output data/interim/deepbgc_predictions.parquet
# ```

# %% [markdown]
# ---
# ## Source 3 — GECCO
#
# GECCO also writes a plain TSV, `<genome>.clusters.tsv` — this is the tool
# whose coordinates actually need converting (see the detour above):

# %%
gecco_df = pd.read_csv(FIXTURES / "gecco_sequence.clusters.tsv", sep="\t")
gecco_df[["sequence_id", "cluster_id", "start", "end", "average_p", "type"]]

# %% [markdown]
# `type` is `"Unknown"` for every row in this real run — GECCO's per-class
# probability columns (`nrp_probability`, `polyketide_probability`, ...) carry
# more signal when that happens, but v1 of the converter keeps `type` as-is
# (matches the schema, no modeling judgment call baked into a parser).
#
# Watch the coordinate conversion happen:

# %%
from convert_gecco_to_parquet import rows_to_regions as gecco_rows_to_regions  # noqa: E402

with open(FIXTURES / "gecco_sequence.clusters.tsv", newline="") as fh:
    gecco_rows = list(csv.DictReader(fh, delimiter="\t"))
gecco_regions = gecco_rows_to_regions(gecco_rows)

first_row, first_region = gecco_rows[0], gecco_regions[0]
print(f"raw TSV      : start={first_row['start']}  end={first_row['end']}")
print(f"converted    : start={first_region.start}  end={first_region.end}")
print("→ start - 1, end unchanged — same conversion pattern as MiBIG in Part 1")

# %% [markdown]
# ### Running it for real
#
# ```bash
# python scripts/convert_gecco_to_parquet.py --inspect <gecco output dir>
# python scripts/convert_gecco_to_parquet.py \
#     --input <gecco output dir> --output data/interim/gecco_predictions.parquet
# ```

# %% [markdown]
# ---
# ## Why this matters: one schema → one scorer
#
# All three converters produce the exact same type. Let's prove it by writing
# each tool's regions through the real parquet writer and reading them back with
# the loader `sharp.evaluate` actually uses:

# %%
from IPython.display import display  # noqa: E402

tool_parquet_paths = {}
for name, regions in [
    ("antiSMASH", antismash_regions),
    ("DeepBGC", deepbgc_regions),
    ("GECCO", gecco_regions),
]:
    tmp_path = Path(f"/tmp/{name.lower()}_predictions_demo.parquet")
    write_predictions_parquet(tmp_path, regions)
    roundtripped = load_predictions_parquet(tmp_path)
    tool_parquet_paths[name] = tmp_path
    print(f"{name:10s} {len(roundtripped)} region(s), "
          f"type={type(roundtripped[0]).__name__}")

# %% [markdown]
# Same schema, side by side — `head()` of each tool's actual
# `predictions.parquet`, read straight off disk with pandas (not through
# `PredictedRegion`, just to eyeball the raw table each `sharp.evaluate` call
# will consume):

# %%
for name, path in tool_parquet_paths.items():
    print(f"--- {name}: {path} ---")
    display(pd.read_parquet(path).head())

# %% [markdown]
# Because every tool produces a list of `PredictedRegion`, the scoring code
# doesn't branch on which tool made the prediction. A tiny sanity check with
# the real metric (`reciprocal_overlap`, from Part 1) makes the point: this
# antiSMASH region, checked against an artificial "known cluster" placed right
# on top of it, matches — and checked against one far away, doesn't — using
# the exact same function that will later score DeepBGC, GECCO, and S(H)ARP
# itself:

# %%
from sharp.io import KnownCluster  # noqa: E402
from sharp.metrics import reciprocal_overlap  # noqa: E402

pred = antismash_regions[0]
print("prediction (antiSMASH):", pred)

on_target = KnownCluster(
    cluster_id="demo", contig=pred.contig, start=pred.start, end=pred.end,
)
off_target = KnownCluster(
    cluster_id="demo", contig=pred.contig,
    start=pred.start + 10_000_000, end=pred.end + 10_000_000,
)

print("on-target known cluster matches? ", reciprocal_overlap(pred, on_target, min_frac=0.5))
print("off-target known cluster matches?", reciprocal_overlap(pred, off_target, min_frac=0.5))

# %% [markdown]
# ---
# ## Recap
#
# - Three baseline tools, three different output formats (JSON, two differently-shaped
#   TSVs) — all normalized into **one 5-column `PredictedRegion` table**.
# - S(H)ARP never runs the tools; each converter only **parses output files** that
#   already exist, in the S(H)ARP env, with no subprocess call.
# - Coordinate conventions were **measured, not assumed** — span vs. matching `.gbk`
#   LOCUS length, checked across every row. antiSMASH and DeepBGC needed no
#   conversion; GECCO needed `start - 1`.
# - Once normalized, `sharp.evaluate` scores **any tool against either ground
#   truth** (Part 1) with the same code — the same design idea that motivated
#   normalizing ground truth in the first place.
#
# **Next:** run `sharp.evaluate` for real — antiSMASH, DeepBGC, GECCO, and S(H)ARP,
# all against the same MiBIG/BGC Atlas ground truth from Part 1 — to get an actual
# benchmark comparison table.
