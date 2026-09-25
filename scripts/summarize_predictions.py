#!/usr/bin/env python3
"""Raw-number comparison of two or more tools over the same genome database.

This is the no-ground-truth counterpart to `sharp.evaluate`. With no MiBiG
labels there is no recall and no precision, only *what each tool called*, so
this reports the quantities that can honestly be compared between tools and
normalises away the two things that otherwise dominate them.

Three warnings that the numbers themselves cannot carry, so read them once:

1. **A region count is not a quality score.** More calls is not better; with no
   labels, nothing here says which tool is right.

2. **Region counts are not commensurable across tools.** antiSMASH merges
   neighbouring protoclusters into one region while DeepBGC emits separate
   candidates, so two tools calling identical territory can differ ~2x in count
   and not at all in `bp_called`. That is why every table reports both, and why
   `frac_bp_called` is the fairer headline.

3. **A count is a threshold, not a fact,** for any tool with a score. DeepBGC's
   `deepbgc_score` and GECCO's `average_p` move the count freely; antiSMASH is
   rule-based and its converter sets `p_bgc = 1.0`, so it is unaffected by
   `--thresholds` and appears once per threshold with identical numbers. Sweep
   the threshold and report the curve, not one point.

Outputs (`--output-dir`):

    summary_by_tool.tsv       one row per (tool, threshold) — the headline
    summary_by_assembly.tsv   one row per (assembly, tool, threshold)

Both are keyed by `assembly` + `contig` via the manifest, which is the join a
collaborator running a third tool needs: give them `genomes.tsv` and ask for
the same two tables.

Usage:
    python scripts/summarize_predictions.py \\
        --manifest data/interim/actino_db/genomes.tsv \\
        --predictions antismash=data/interim/antismash_predictions_actino.parquet \\
                      deepbgc=data/interim/deepbgc_predictions_actino.parquet \\
        --thresholds 0.0 0.5 0.8 \\
        --output-dir data/processed/actino_comparison
"""
from __future__ import annotations

import argparse
import logging
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from build_genome_manifest import load_contig_index  # noqa: E402
from sharp.io import PredictedRegion, load_predictions_parquet  # noqa: E402
from sharp.metrics import merge_intervals  # noqa: E402

LOG = logging.getLogger("summarize_predictions")


@dataclass(frozen=True)
class ToolSummary:
    """One tool at one threshold, over the whole database."""
    tool: str
    threshold: float
    n_regions: int
    n_assemblies_called: int
    bp_called: int
    total_bp: int
    median_region_bp: int

    @property
    def regions_per_mb(self) -> float:
        return self.n_regions / (self.total_bp / 1e6) if self.total_bp else 0.0

    @property
    def frac_bp_called(self) -> float:
        return self.bp_called / self.total_bp if self.total_bp else 0.0


@dataclass(frozen=True)
class AssemblySummary:
    """One tool at one threshold, over one assembly."""
    assembly: str
    tool: str
    threshold: float
    n_regions: int
    bp_called: int
    genome_bp: int

    @property
    def regions_per_mb(self) -> float:
        return self.n_regions / (self.genome_bp / 1e6) if self.genome_bp else 0.0

    @property
    def frac_bp_called(self) -> float:
        return self.bp_called / self.genome_bp if self.genome_bp else 0.0


def parse_spec(spec: str) -> tuple[str, Path]:
    """`antismash=path/to.parquet` -> ("antismash", Path(...))."""
    if "=" not in spec:
        raise SystemExit(
            f"--predictions takes NAME=PATH, got {spec!r} "
            "(e.g. antismash=data/interim/antismash_predictions.parquet)"
        )
    name, _, path = spec.partition("=")
    if not name or not path:
        raise SystemExit(f"--predictions takes NAME=PATH, got {spec!r}")
    return name, Path(path)


def bp_called_per_contig(regions: list[PredictedRegion]) -> dict[str, int]:
    """Merged bp called on each contig.

    Merging first is what makes this comparable between tools: overlapping or
    adjacent calls from a splitter are counted once, the same as one call from
    a merger covering the same territory.
    """
    by_contig: dict[str, list[tuple[int, int]]] = {}
    for r in regions:
        by_contig.setdefault(r.contig, []).append((r.start, r.end))
    return {
        contig: sum(e - s for s, e in merge_intervals(intervals))
        for contig, intervals in by_contig.items()
    }


def summarize(
    tool: str,
    regions: list[PredictedRegion],
    index: dict[str, tuple[str, int]],
    threshold: float,
) -> tuple[ToolSummary, list[AssemblySummary]]:
    """Summarise one tool's predictions at one score threshold."""
    kept = [r for r in regions if r.p_bgc >= threshold]

    # Total bp is the whole database, not just the contigs that got a call:
    # a genome a tool found nothing in still counts against its calls-per-Mb.
    total_bp = sum(length for _, length in index.values())
    genome_bp: dict[str, int] = {}
    for assembly, length in index.values():
        genome_bp[assembly] = genome_bp.get(assembly, 0) + length

    per_contig_bp = bp_called_per_contig(kept)
    n_by_assembly: dict[str, int] = {}
    bp_by_assembly: dict[str, int] = {}
    for r in kept:
        assembly = index[r.contig][0]
        n_by_assembly[assembly] = n_by_assembly.get(assembly, 0) + 1
    for contig, bp in per_contig_bp.items():
        assembly = index[contig][0]
        bp_by_assembly[assembly] = bp_by_assembly.get(assembly, 0) + bp

    rows = [
        AssemblySummary(
            assembly=assembly,
            tool=tool,
            threshold=threshold,
            n_regions=n_by_assembly.get(assembly, 0),
            bp_called=bp_by_assembly.get(assembly, 0),
            genome_bp=bp,
        )
        for assembly, bp in sorted(genome_bp.items())
    ]

    spans = sorted(r.end - r.start for r in kept)
    overall = ToolSummary(
        tool=tool,
        threshold=threshold,
        n_regions=len(kept),
        n_assemblies_called=sum(1 for r in rows if r.n_regions > 0),
        bp_called=sum(per_contig_bp.values()),
        total_bp=total_bp,
        median_region_bp=int(statistics.median(spans)) if spans else 0,
    )
    return overall, rows


def check_contigs(tool: str, regions: list[PredictedRegion],
                  index: dict[str, tuple[str, int]]) -> list[PredictedRegion]:
    """Drop regions on contigs the manifest does not know, loudly.

    This is the failure mode worth catching early: if a tool renames or
    truncates contig identifiers (dropping the `.1` version suffix, say),
    every region silently fails to join and the tool looks like it found
    nothing. A large `unknown` count means the names do not match, not that
    the tool was quiet.
    """
    known = [r for r in regions if r.contig in index]
    unknown = [r for r in regions if r.contig not in index]
    if unknown:
        names = sorted({r.contig for r in unknown})
        LOG.warning(
            "%s: %d region(s) on %d contig(s) absent from the manifest — "
            "dropped. First few: %s",
            tool, len(unknown), len(names), ", ".join(names[:5]),
        )
        if len(known) == 0:
            LOG.error(
                "%s: NOTHING joined to the manifest. The contig identifiers "
                "almost certainly disagree (version suffix? assembly name "
                "instead of accession?) — fix that before reading any number.",
                tool,
            )
    return known


# ══════════════════════════════ writing ════════════════════════════════════

TOOL_HEADER = ("tool\tthreshold\tn_regions\tn_assemblies_called\tbp_called\t"
               "total_bp\tregions_per_Mb\tfrac_bp_called\tmedian_region_bp\n")
ASM_HEADER = ("assembly\ttool\tthreshold\tn_regions\tbp_called\tgenome_bp\t"
              "regions_per_Mb\tfrac_bp_called\n")


def write_tool_tsv(path: Path, rows: list[ToolSummary]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(TOOL_HEADER)
        for r in rows:
            fh.write(f"{r.tool}\t{r.threshold:g}\t{r.n_regions}\t"
                     f"{r.n_assemblies_called}\t{r.bp_called}\t{r.total_bp}\t"
                     f"{r.regions_per_mb:.3f}\t{r.frac_bp_called:.5f}\t"
                     f"{r.median_region_bp}\n")
    return len(rows)


def write_assembly_tsv(path: Path, rows: list[AssemblySummary]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write(ASM_HEADER)
        for r in rows:
            fh.write(f"{r.assembly}\t{r.tool}\t{r.threshold:g}\t{r.n_regions}\t"
                     f"{r.bp_called}\t{r.genome_bp}\t{r.regions_per_mb:.3f}\t"
                     f"{r.frac_bp_called:.5f}\n")
    return len(rows)


def print_table(rows: list[ToolSummary]) -> None:
    print(f"\n{'='*86}\nRAW COMPARISON — no ground truth, so no recall and no "
          f"precision\n{'='*86}")
    print(f"{'tool':<14} {'thr':>5} {'regions':>9} {'genomes':>8} "
          f"{'reg/Mb':>8} {'bp called':>14} {'%genome':>8} {'median bp':>10}")
    for r in rows:
        print(f"{r.tool:<14} {r.threshold:>5.2f} {r.n_regions:>9,} "
              f"{r.n_assemblies_called:>8,} {r.regions_per_mb:>8.3f} "
              f"{r.bp_called:>14,} {100*r.frac_bp_called:>7.2f}% "
              f"{r.median_region_bp:>10,}")
    print()


# ══════════════════════════════ orchestration ══════════════════════════════

def run(
    manifest_path: Path,
    specs: list[str],
    thresholds: list[float],
    output_dir: Path | None,
) -> None:
    index = load_contig_index(manifest_path)
    LOG.info("manifest: %d contigs across %d assemblies, %.1f Gb",
             len(index), len({a for a, _ in index.values()}),
             sum(n for _, n in index.values()) / 1e9)

    tool_rows: list[ToolSummary] = []
    asm_rows: list[AssemblySummary] = []

    for spec in specs:
        tool, path = parse_spec(spec)
        if not path.is_file():
            raise SystemExit(f"no such predictions file: {path}")
        regions = check_contigs(tool, load_predictions_parquet(path), index)
        LOG.info("%s: %d region(s) joined to the manifest", tool, len(regions))
        for t in thresholds:
            overall, rows = summarize(tool, regions, index, t)
            tool_rows.append(overall)
            asm_rows.extend(rows)

    print_table(tool_rows)

    if output_dir is None:
        LOG.info("no --output-dir given; nothing written")
        return
    n_t = write_tool_tsv(output_dir / "summary_by_tool.tsv", tool_rows)
    n_a = write_assembly_tsv(output_dir / "summary_by_assembly.tsv", asm_rows)
    LOG.info("wrote %d tool row(s) and %d assembly row(s) → %s",
             n_t, n_a, output_dir)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--manifest", type=Path, required=True,
                   help="genomes.tsv from build_genome_manifest.py")
    p.add_argument("--predictions", nargs="+", required=True, metavar="NAME=PATH",
                   help="one or more labelled predictions parquets, "
                        "e.g. antismash=data/interim/antismash_predictions.parquet")
    p.add_argument("--thresholds", nargs="+", type=float, default=[0.0],
                   metavar="P",
                   help="p_bgc cutoffs to report (default: 0.0, i.e. keep "
                        "everything). antiSMASH is rule-based and unaffected.")
    p.add_argument("--output-dir", type=Path, default=None,
                   help="write summary_by_tool.tsv and summary_by_assembly.tsv here")
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
        manifest_path=args.manifest,
        specs=args.predictions,
        thresholds=sorted(args.thresholds),
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
