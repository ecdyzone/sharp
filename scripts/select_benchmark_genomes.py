#!/usr/bin/env python3
"""Select the genome set for a scaled-up baseline benchmark.

Turns a ground-truth TSV into three files that together define one benchmark
run: which genomes to fetch, which contigs are in scope, and the ground truth
restricted to that scope.

Why this is not just "sort by cluster count":

  1. Most ground-truth records are not genomes. ~58% of the MiBIG 4.0
     Streptomyces records are BGC-only deposits — the record *is* the cluster
     (BGC0001709 spans 0..119,206 of a 119 kb record). A detector run on one
     of those scores detection recall ~1.0 by construction, for every tool, so
     mixing them into a headline table flatters all baselines equally and
     measures nothing. `--min-length` drops them.

  2. The same physical sequence appears under several accessions. RefSeq and
     GenBank copies are byte-identical but carry different accessions, and
     MiBIG files clusters against whichever the submitter used. NC_003888.3 and
     AL645882.2 are both the 8,667,507 bp S. coelicolor A3(2) chromosome, with
     15 clusters filed under one and 1 under the other. Run naively you either
     analyze the genome twice or lose that 16th cluster from the denominator.

Point 2 is why this script emits a ground truth of its own. Merged records are
collapsed onto a single primary accession, and the ground truth is rewritten to
match — coordinates are unchanged, because the merged accessions denote the
same sequence. Emitting a scope file without rewriting the ground truth would
silently drop every cluster filed under a non-primary accession.

Record lengths come from NCBI esummary and are cached, so a rerun (and the test
suite) needs no network. Only lengths are fetched — not sequences.

Usage:
    # Inspect the selection: ranking, what was filtered, what was merged
    python scripts/select_benchmark_genomes.py \\
        --ground-truth data/raw/streptomyces_ground_truth.tsv --inspect

    # Write the benchmark set (default: top 50 genome-scale records)
    python scripts/select_benchmark_genomes.py \\
        --ground-truth data/raw/streptomyces_ground_truth.tsv \\
        --output-dir data/interim/benchmark_set

    # Take every genome-scale record instead of the top N
    python scripts/select_benchmark_genomes.py \\
        --ground-truth data/raw/streptomyces_ground_truth.tsv \\
        --output-dir data/interim/benchmark_set --top-n 0

Outputs (into --output-dir):
    benchmark_genomes.tsv       one row per selected genome, with provenance
    analyzed_contigs.txt        the --contigs scope file for sharp.evaluate
    benchmark_ground_truth.tsv  ground truth restricted to the selection, with
                                contig names normalized to primary accessions
"""
from __future__ import annotations

import argparse
import json
import logging
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from sharp.config import INTERIM_DIR, RAW_DIR
from sharp.io import KnownCluster, load_ground_truth_tsv, write_ground_truth_tsv

LOG = logging.getLogger("select_benchmark_genomes")

#: Records shorter than this are BGC-only deposits or short contigs, not
#: genomes. 1 Mb sits well above the largest MiBIG locus (~150 kb) and below
#: the smallest real Streptomyces replicon in the set (the 1.37 Mb
#: S. ambofaciens chromosomal arm and the 1.79 Mb S. clavuligerus plasmid
#: pSCL4, both of which we want to keep).
DEFAULT_MIN_LENGTH = 1_000_000

DEFAULT_TOP_N = 50

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


@dataclass(frozen=True)
class RecordInfo:
    """What NCBI knows about one nucleotide accession."""

    accession: str  # versioned, as NCBI reports it
    length: int
    title: str
    #: False when `length` is a lower bound derived from ground-truth
    #: coordinates rather than reported by NCBI. See `infer_lower_bounds`.
    length_known: bool = True


@dataclass(frozen=True)
class BenchmarkGenome:
    """One genome to run the baselines on, after merging equivalent records."""

    accession: str  # the primary accession — what we fetch and analyze
    accessions: tuple[str, ...]  # every accession merged into this genome
    length: int
    organism: str
    clusters: tuple[KnownCluster, ...]
    length_known: bool = True

    @property
    def n_clusters(self) -> int:
        return len(self.clusters)


# ═══════════════════════════ record length lookup ══════════════════════════
# The only network-touching part. Isolated here so every rule below stays pure
# and testable offline.
# ═══════════════════════════════════════════════════════════════════════════

def read_length_cache(path: Path) -> dict[str, RecordInfo]:
    """Load a previously written length cache. Missing file → empty dict."""
    if not path.exists():
        return {}
    out: dict[str, RecordInfo] = {}
    with path.open() as fh:
        header = fh.readline()
        if not header.startswith("caption"):
            fh.seek(0)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            # caption, accession, length, title — caption is the lookup key and
            # is deliberately separate from the versioned accession, because
            # ground truth and NCBI disagree on version suffixes.
            if len(parts) < 4 or not parts[2].isdigit():
                continue
            out[parts[0]] = RecordInfo(parts[1], int(parts[2]), parts[3])
    return out


def write_length_cache(path: Path, info: dict[str, RecordInfo]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("caption\taccession\tlength\ttitle\n")
        for caption in sorted(info):
            r = info[caption]
            fh.write(f"{caption}\t{r.accession}\t{r.length}\t{r.title}\n")


def fetch_record_lengths(
    accessions: list[str], batch_size: int = 120, pause: float = 0.4
) -> dict[str, RecordInfo]:
    """Fetch `slen` and `title` for each accession via NCBI esummary.

    esummary accepts accessions directly as `id`, which is why this does not
    esearch first — and notably why it works for accessions the `[Accession]`
    search field fails to resolve (BA000030.3 among them).

    Keyed by the *unversioned* caption, because ground truth and NCBI do not
    always agree on the version suffix (the ground truth has bare `CP114200`
    where NCBI reports `CP114200.1`).
    """
    out: dict[str, RecordInfo] = {}
    for i in range(0, len(accessions), batch_size):
        batch = accessions[i : i + batch_size]
        params = urllib.parse.urlencode(
            {"db": "nuccore", "retmode": "json", "id": ",".join(batch)}
        )
        LOG.info("esummary batch %d–%d of %d", i + 1, i + len(batch), len(accessions))
        try:
            with urllib.request.urlopen(EUTILS + "esummary.fcgi?" + params, timeout=60) as fh:
                data = json.load(fh)
        except Exception as e:  # noqa: BLE001 — network is best-effort here
            LOG.warning("esummary batch failed (%s); those records stay unresolved", e)
            continue
        result = data.get("result", {})
        for uid in result.get("uids", []):
            rec = result[uid]
            slen = rec.get("slen")
            if not isinstance(slen, int):
                continue
            caption = rec.get("caption") or ""
            out[caption] = RecordInfo(
                accession=rec.get("accessionversion") or caption,
                length=slen,
                title=rec.get("title", ""),
            )
        time.sleep(pause)
    return out


def resolve_lengths(
    accessions: list[str], cache_path: Path, offline: bool = False
) -> dict[str, RecordInfo]:
    """Return caption → RecordInfo, fetching only what the cache lacks."""
    cache = read_length_cache(cache_path)
    captions = {a: caption_of(a) for a in accessions}
    missing = sorted({c for c in captions.values() if c not in cache})
    if missing and not offline:
        LOG.info("%d/%d accessions not cached — querying NCBI",
                 len(missing), len(set(captions.values())))
        fetched = fetch_record_lengths(missing)
        cache.update(fetched)
        write_length_cache(cache_path, cache)
        LOG.info("length cache now holds %d records → %s", len(cache), cache_path)
    elif missing:
        LOG.warning("--offline: %d accessions have no cached length and will "
                    "be dropped", len(missing))
    return cache


def infer_lower_bounds(
    clusters: list[KnownCluster], info: dict[str, RecordInfo]
) -> dict[str, RecordInfo]:
    """Length lower bounds for accessions NCBI could not resolve.

    Some WGS contigs are absent from the nuccore esearch index and are rejected
    by esummary as "Invalid uid" (JAFMOF010000002.1, JAAVVL010000001.1 and
    friends). Rather than drop them silently, fall back on the ground truth: a
    cluster ending at 1,725,178 proves its record is at least that long.

    The bound can only under-estimate, so it never promotes a BGC-only deposit
    into the genome-scale set — it only risks excluding a genuine genome whose
    single known cluster sits near the start, which is the safe direction to
    err. Records resolved this way are never merged with a twin, since twin
    detection needs an exact length.
    """
    bounds: dict[str, int] = {}
    for c in clusters:
        cap = caption_of(c.contig)
        if cap in info:
            continue
        bounds[cap] = max(bounds.get(cap, 0), c.end)
    return {
        cap: RecordInfo(accession=cap, length=end, title=cap, length_known=False)
        for cap, end in bounds.items()
    }


# ═══════════════════════════ selection rules (pure) ════════════════════════

def caption_of(accession: str) -> str:
    """Strip the version suffix. `CP114200.1` → `CP114200`, `CP114200` → same."""
    return accession.split(".")[0]


def organism_key(title: str) -> str:
    """Collapse an NCBI title to the strain that owns it.

    The first three words carry genus, species, and strain — enough to tell
    `Streptomyces sp. CS113` from `Streptomyces sp. CS147`, while ignoring the
    trailing description that differs between a RefSeq and GenBank copy of one
    sequence ("... genomic scaffold scaffold00001" vs "... scaffold00001").
    """
    return " ".join(title.split()[:3]).rstrip(",").lower()


def group_into_genomes(
    clusters: list[KnownCluster], info: dict[str, RecordInfo]
) -> tuple[list[BenchmarkGenome], dict[str, list[str]]]:
    """Group ground-truth clusters into physical genomes.

    Two accessions denote the same genome when they agree on both length and
    organism — the test that identifies RefSeq/GenBank twins (NC_003888.3 and
    AL645882.2, both 8,667,507 bp of S. coelicolor A3(2)) without relying on
    any textual relationship between the accessions, of which there is none.

    Returns (genomes, unresolved) where `unresolved` maps a rejection reason to
    the contigs it dropped.
    """
    unresolved: dict[str, list[str]] = {}
    buckets: dict[tuple[int, str], list[KnownCluster]] = {}
    bucket_accs: dict[tuple[int, str], set[str]] = {}
    bucket_info: dict[tuple[int, str], RecordInfo] = {}

    for c in clusters:
        rec = info.get(caption_of(c.contig))
        if rec is None:
            unresolved.setdefault("no length from NCBI", []).append(c.contig)
            continue
        # An exact length is what makes twin detection safe; a lower-bound
        # record gets a bucket of its own so it can never be merged.
        key = ((rec.length, organism_key(rec.title)) if rec.length_known
               else (-1, rec.accession))
        buckets.setdefault(key, []).append(c)
        bucket_accs.setdefault(key, set()).add(rec.accession)
        bucket_info[key] = rec

    genomes: list[BenchmarkGenome] = []
    for key, rows in buckets.items():
        rec = bucket_info[key]
        accs = tuple(sorted(bucket_accs[key]))
        genomes.append(BenchmarkGenome(
            accession=primary_accession(accs),
            accessions=accs,
            length=rec.length,
            organism=" ".join(rec.title.split()[:3]).rstrip(","),
            clusters=tuple(rows),
            length_known=rec.length_known,
        ))
    return genomes, unresolved


def primary_accession(accessions: tuple[str, ...]) -> str:
    """Pick which accession of a merged group to fetch and analyze.

    Prefer an INSDC accession over a RefSeq one (`NC_`/`NZ_`): ground truth is
    overwhelmingly INSDC-keyed, so choosing it keeps the rewrite in
    `normalize_ground_truth` to a minimum. Ties break alphabetically so the
    choice is deterministic across runs.
    """
    insdc = [a for a in accessions if not a.upper().startswith(("NC_", "NZ_"))]
    return sorted(insdc or list(accessions))[0]


def select_genomes(
    genomes: list[BenchmarkGenome],
    top_n: int = DEFAULT_TOP_N,
    min_length: int = DEFAULT_MIN_LENGTH,
) -> tuple[list[BenchmarkGenome], list[BenchmarkGenome]]:
    """Filter to genome-scale records and rank. Returns (selected, rejected).

    Ranked by cluster count first — clusters are the scarce resource, and a
    record carrying 15 of them is worth far more than one carrying 1 — then by
    length, which breaks the large tie among single-cluster records in favour
    of complete genomes over short scaffolds.
    """
    kept = [g for g in genomes if g.length >= min_length]
    rejected = [g for g in genomes if g.length < min_length]
    kept.sort(key=lambda g: (-g.n_clusters, -g.length, g.accession))
    if top_n > 0:
        kept = kept[:top_n]
    return kept, rejected


def normalize_ground_truth(selected: list[BenchmarkGenome]) -> list[KnownCluster]:
    """Ground truth for the selection, with contigs renamed to the primary
    accession of the genome they sit on.

    Coordinates are untouched: the merged accessions denote the same sequence,
    so an interval valid on one is valid on the other.
    """
    out: list[KnownCluster] = []
    for g in selected:
        for c in g.clusters:
            out.append(KnownCluster(
                cluster_id=c.cluster_id,
                contig=g.accession,
                start=c.start,
                end=c.end,
                cluster_class=c.cluster_class,
            ))
    out.sort(key=lambda c: (c.contig, c.start, c.cluster_id))
    return out


# ══════════════════════════════ output ═════════════════════════════════════

def write_benchmark_set(path: Path, selected: list[BenchmarkGenome]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("rank\taccession\tall_accessions\tlength\tlength_known\t"
                 "n_clusters\tcluster_ids\torganism\n")
        for i, g in enumerate(selected, 1):
            ids = ",".join(c.cluster_id for c in g.clusters)
            fh.write(f"{i}\t{g.accession}\t{','.join(g.accessions)}\t{g.length}\t"
                     f"{str(g.length_known).lower()}\t{g.n_clusters}\t{ids}\t"
                     f"{g.organism}\n")


def write_contigs_file(path: Path, selected: list[BenchmarkGenome]) -> None:
    """The `--contigs` scope file for sharp.evaluate: one primary accession per
    line, matching both the FASTA we fetch and the normalized ground truth."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{g.accession}\n" for g in selected))


# ══════════════════════════════ inspect mode ═══════════════════════════════

def inspect(
    genomes: list[BenchmarkGenome],
    selected: list[BenchmarkGenome],
    rejected: list[BenchmarkGenome],
    unresolved: dict[str, list[str]],
) -> None:
    total_clusters = sum(g.n_clusters for g in genomes)
    print(f"\n{'='*78}\nSELECTION SUMMARY\n{'='*78}")
    print(f"physical genomes after merging equivalent accessions : {len(genomes)}")
    print(f"  clusters on them                                   : {total_clusters}")
    print(f"below --min-length (BGC-only deposits, short contigs) : {len(rejected)} "
          f"records / {sum(g.n_clusters for g in rejected)} clusters")
    print(f"selected                                             : {len(selected)} "
          f"records / {sum(g.n_clusters for g in selected)} clusters")
    for reason, contigs in sorted(unresolved.items()):
        print(f"unresolved ({reason}): {len(contigs)} — {', '.join(sorted(set(contigs))[:6])}")

    approx = [g for g in selected if not g.length_known]
    if approx:
        print(f"\n{'-'*78}\nLENGTH IS A LOWER BOUND (NCBI could not resolve these)\n{'-'*78}")
        for g in approx:
            print(f"  {g.accession:<20} >= {g.length:>11,} bp  "
                  f"{g.n_clusters} cluster(s) — not eligible for twin merging")

    merged = [g for g in genomes if len(g.accessions) > 1]
    if merged:
        print(f"\n{'-'*78}\nMERGED RECORDS (same length + organism = same sequence)\n{'-'*78}")
        for g in sorted(merged, key=lambda x: -x.n_clusters):
            print(f"  {g.length:>11,} bp  {g.n_clusters:>2} clusters  "
                  f"{' = '.join(g.accessions)}  [{g.organism}]")

    print(f"\n{'-'*78}\nSELECTED GENOMES\n{'-'*78}")
    print(f"{'rank':>4} {'clus':>4} {'length':>12}  {'accession':<16} organism")
    for i, g in enumerate(selected, 1):
        print(f"{i:>4} {g.n_clusters:>4} {g.length:>12,}  {g.accession:<16} {g.organism}")
    print(f"\ntotal sequence to analyze: {sum(g.length for g in selected):,} bp")


# ══════════════════════════════ orchestration ══════════════════════════════

def run(
    ground_truth: Path,
    output_dir: Path | None,
    top_n: int,
    min_length: int,
    cache_path: Path,
    offline: bool,
    do_inspect: bool,
) -> None:
    clusters = load_ground_truth_tsv(ground_truth)
    LOG.info("loaded %d ground-truth clusters from %s", len(clusters), ground_truth)

    accessions = sorted({c.contig for c in clusters})
    info = resolve_lengths(accessions, cache_path, offline=offline)

    bounds = infer_lower_bounds(clusters, info)
    if bounds:
        LOG.warning("%d accession(s) unresolvable at NCBI — using a length "
                    "lower bound from ground-truth coordinates: %s",
                    len(bounds), ", ".join(sorted(bounds)[:6])
                    + (" ..." if len(bounds) > 6 else ""))
        info = {**info, **bounds}

    genomes, unresolved = group_into_genomes(clusters, info)
    selected, rejected = select_genomes(genomes, top_n=top_n, min_length=min_length)

    for reason, contigs in sorted(unresolved.items()):
        LOG.warning("%d contig(s) unresolved — %s", len(set(contigs)), reason)
    LOG.info("%d physical genomes; %d below --min-length; selected %d "
             "carrying %d clusters",
             len(genomes), len(rejected), len(selected),
             sum(g.n_clusters for g in selected))

    if do_inspect:
        inspect(genomes, selected, rejected, unresolved)
        return

    if output_dir is None:
        raise SystemExit("--output-dir is required unless --inspect is given")
    output_dir.mkdir(parents=True, exist_ok=True)
    write_benchmark_set(output_dir / "benchmark_genomes.tsv", selected)
    write_contigs_file(output_dir / "analyzed_contigs.txt", selected)
    gt = normalize_ground_truth(selected)
    write_ground_truth_tsv(output_dir / "benchmark_ground_truth.tsv", gt)

    LOG.info("wrote %d genomes → %s", len(selected),
             output_dir / "benchmark_genomes.tsv")
    LOG.info("wrote %d contigs → %s", len(selected),
             output_dir / "analyzed_contigs.txt")
    LOG.info("wrote %d clusters → %s", len(gt),
             output_dir / "benchmark_ground_truth.tsv")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--ground-truth", type=Path,
                   default=RAW_DIR / "streptomyces_ground_truth.tsv",
                   help="ground-truth TSV to select from (default: %(default)s)")
    p.add_argument("--output-dir", type=Path, default=None,
                   help="directory for the three output files "
                        "(required unless --inspect)")
    p.add_argument("--top-n", type=int, default=DEFAULT_TOP_N,
                   help="keep the N highest-ranked genomes; 0 keeps all "
                        "genome-scale records (default: %(default)s)")
    p.add_argument("--min-length", type=int, default=DEFAULT_MIN_LENGTH,
                   help="drop records shorter than this many bp — they are "
                        "BGC-only deposits, not genomes (default: %(default)s)")
    p.add_argument("--lengths-cache", type=Path,
                   default=INTERIM_DIR / "record_lengths.tsv",
                   help="cache of NCBI record lengths (default: %(default)s)")
    p.add_argument("--offline", action="store_true",
                   help="never query NCBI; drop accessions missing from the cache")
    p.add_argument("--inspect", action="store_true",
                   help="print the ranking, merges and filtering, and write nothing")
    p.add_argument("-v", "--verbose", action="store_true")
    return p


def main() -> None:
    args = build_parser().parse_args()
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )
    run(
        ground_truth=args.ground_truth,
        output_dir=args.output_dir,
        top_n=args.top_n,
        min_length=args.min_length,
        cache_path=args.lengths_cache,
        offline=args.offline,
        do_inspect=args.inspect,
    )


if __name__ == "__main__":
    main()
