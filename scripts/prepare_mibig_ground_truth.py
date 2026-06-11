#!/usr/bin/env python3
"""Build a ground-truth TSV from a MiBIG JSON dump.

Walks a directory of MiBIG per-entry JSON files (e.g. the contents of
mibig_json_4.0.tar.gz) and emits a TSV with one row per cluster locus,
in the format `sharp.io.load_ground_truth_tsv` expects:

    cluster_id  contig  start  end  class

IMPORTANT — verify the schema first:
    The MiBIG JSON schema changed across major versions. Before trusting
    the output, run:

        python scripts/prepare_mibig_ground_truth.py --inspect data/raw/mibig_json_4.0

    This prints the structure of a few real entries so you can confirm the
    field paths this script reads. If they differ, the FIELD PATHS section
    below is the only thing that needs editing.

Usage:
    # Inspect the real schema (do this first)
    python scripts/prepare_mibig_ground_truth.py --inspect data/raw/mibig_json_4.0

    # Build the TSV
    python scripts/prepare_mibig_ground_truth.py \\
        --input-dir data/raw/mibig_json_4.0 \\
        --output data/raw/mibig_ground_truth.tsv

    # Restrict to a genus (e.g. only Streptomyces) for a focused benchmark
    python scripts/prepare_mibig_ground_truth.py \\
        --input-dir data/raw/mibig_json_4.0 \\
        --output data/raw/streptomyces_ground_truth.tsv \\
        --genus Streptomyces
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Any

from sharp.config import RAW_DIR
from sharp.io import KnownCluster, write_ground_truth_tsv

LOG = logging.getLogger("prepare_mibig_ground_truth")


# ══════════════════════════════ FIELD PATHS ════════════════════════════════
# These helpers isolate every assumption about the MiBIG JSON layout in ONE
# place. If `--inspect` shows different paths, edit only this section.
#
# As of MiBIG 4.0, a typical entry looks like (abbreviated):
#   {
#     "accession": "BGC0000001",            # or nested under "cluster"
#     "biosynthesis": {"classes": [{"class": "PKS"}, ...]},
#     "loci": [
#        {"accession": "AB000001.1",
#         "location": {"from": 1, "to": 41001}}
#     ],
#     "taxonomy": {"name": "Streptomyces ..."}     # genus lives in the name
#   }
# Older (3.x) entries nest everything under a top-level "cluster" key and use
# "mibig_accession", "biosyn_class", and "loci.accession/start/end". The
# accessors below try the 4.0 path first, then fall back to 3.x.
# ════════════════════════════════════════════════════════════════════════════

def _root(entry: dict[str, Any]) -> dict[str, Any]:
    """3.x wraps everything under 'cluster'; 4.0 is flat. Return the root
    that holds the cluster fields either way."""
    return entry.get("cluster", entry)


def get_cluster_id(entry: dict[str, Any], fallback: str) -> str:
    r = _root(entry)
    return (
        r.get("accession")
        or r.get("mibig_accession")
        or entry.get("accession")
        or fallback
    )


def get_classes(entry: dict[str, Any]) -> list[str]:
    r = _root(entry)
    # 4.0: biosynthesis.classes -> [{"class": "PKS"}, ...]
    bio = r.get("biosynthesis")
    if isinstance(bio, dict):
        classes = bio.get("classes")
        if isinstance(classes, list):
            out = []
            for c in classes:
                if isinstance(c, dict) and "class" in c:
                    out.append(str(c["class"]))
                elif isinstance(c, str):
                    out.append(c)
            if out:
                return out
    # 3.x: biosyn_class -> ["Polyketide", ...]
    legacy = r.get("biosyn_class")
    if isinstance(legacy, list):
        return [str(x) for x in legacy]
    return []


def get_taxonomy_name(entry: dict[str, Any]) -> str | None:
    r = _root(entry)
    tax = r.get("taxonomy")
    if isinstance(tax, dict):
        return tax.get("name") or tax.get("organism")
    # 3.x sometimes: organism_name at root
    return r.get("organism_name") or r.get("organism")


def get_loci(entry: dict[str, Any]) -> list[dict[str, Any]]:
    """Return a list of locus dicts. 4.0 uses a 'loci' list; some 3.x entries
    use a single 'loci' object."""
    r = _root(entry)
    loci = r.get("loci")
    if isinstance(loci, list):
        return loci
    if isinstance(loci, dict):
        return [loci]
    return []


def get_locus_coords(locus: dict[str, Any]) -> tuple[str | None, int | None, int | None]:
    """Extract (accession/contig, start, end) from one locus dict.

    Returns (None, None, None) components that are missing — caller decides
    how to handle partial loci.
    """
    contig = locus.get("accession") or locus.get("genbank_accession")

    # 4.0: location.from / location.to
    loc = locus.get("location")
    if isinstance(loc, dict):
        start = loc.get("from", loc.get("start"))
        end = loc.get("to", loc.get("end"))
    else:
        # 3.x: start/end directly on the locus
        start = locus.get("start_coord", locus.get("start"))
        end = locus.get("end_coord", locus.get("end"))

    start = int(start) if start is not None else None
    end = int(end) if end is not None else None
    return contig, start, end


# ══════════════════════════════ parsing ════════════════════════════════════

def entry_to_clusters(
    entry: dict[str, Any], source_filename: str
) -> list[KnownCluster]:
    """Convert one MiBIG JSON entry into KnownCluster rows (one per locus
    with usable coordinates). Returns [] if no locus has coordinates."""
    cluster_id = get_cluster_id(entry, fallback=Path(source_filename).stem)
    classes = get_classes(entry)
    cluster_class = "/".join(classes) if classes else None

    out: list[KnownCluster] = []
    loci = get_loci(entry)
    multi = len(loci) > 1
    for i, locus in enumerate(loci):
        contig, start, end = get_locus_coords(locus)
        if contig is None or start is None or end is None:
            continue
        # Half-open convention: MiBIG coords are 1-based inclusive [from, to].
        # Convert to 0-based half-open [start-1, end) used everywhere else.
        norm_start = max(0, start - 1)
        norm_end = end
        if norm_end <= norm_start:
            continue
        # If an entry has multiple loci, suffix the id to keep ids unique.
        cid = f"{cluster_id}.{i+1}" if multi else cluster_id
        out.append(KnownCluster(
            cluster_id=cid,
            contig=contig,
            start=norm_start,
            end=norm_end,
            cluster_class=cluster_class,
        ))
    return out


def load_entries(input_dir: Path) -> list[tuple[dict[str, Any], str]]:
    """Read every .json file in input_dir (recursively). Returns a list of
    (parsed_json, filename) tuples. Malformed files are skipped with a warning."""
    files = sorted(input_dir.rglob("*.json"))
    out: list[tuple[dict[str, Any], str]] = []
    for f in files:
        try:
            out.append((json.loads(f.read_text()), f.name))
        except (json.JSONDecodeError, OSError) as e:
            LOG.warning("skip %s — %s", f.name, e)
    return out


# ══════════════════════════════ inspect mode ═══════════════════════════════

def inspect(input_dir: Path, n: int = 3) -> None:
    """Print the structure of the first n entries so the user can verify the
    field paths this script reads against the actual schema on disk."""
    files = sorted(input_dir.rglob("*.json"))
    if not files:
        LOG.error("no .json files found under %s", input_dir)
        return

    LOG.info("found %d JSON files; inspecting first %d", len(files), min(n, len(files)))
    for f in files[:n]:
        entry = json.loads(f.read_text())
        print(f"\n{'='*70}\nFILE: {f.name}")
        print(f"{'='*70}")
        print("top-level keys:", sorted(entry.keys()))
        r = _root(entry)
        if r is not entry:
            print("(everything nested under 'cluster')")
            print("cluster keys:", sorted(r.keys()))

        # Show what each accessor resolves to — this is the verification.
        print("\n→ accessors resolve to:")
        print(f"   cluster_id    = {get_cluster_id(entry, f.stem)!r}")
        print(f"   classes       = {get_classes(entry)!r}")
        print(f"   taxonomy_name = {get_taxonomy_name(entry)!r}")
        loci = get_loci(entry)
        print(f"   n_loci        = {len(loci)}")
        for j, locus in enumerate(loci[:3]):
            print(f"   locus[{j}] keys = {sorted(locus.keys())}")
            print(f"   locus[{j}] -> coords = {get_locus_coords(locus)}")

        # Show the resulting rows.
        clusters = entry_to_clusters(entry, f.name)
        print(f"\n→ would produce {len(clusters)} row(s):")
        for c in clusters[:3]:
            print(f"   {c}")


# ══════════════════════════════ orchestration ══════════════════════════════

def build_ground_truth(
    input_dir: Path,
    output_path: Path,
    genus: str | None = None,
) -> int:
    entries = load_entries(input_dir)
    if not entries:
        LOG.error("no readable JSON files under %s", input_dir)
        return 0
    LOG.info("loaded %d entries", len(entries))

    clusters: list[KnownCluster] = []
    n_no_coords = 0
    n_filtered_genus = 0
    for entry, fname in entries:
        if genus is not None:
            tax = get_taxonomy_name(entry) or ""
            if genus.lower() not in tax.lower():
                n_filtered_genus += 1
                continue
        rows = entry_to_clusters(entry, fname)
        if not rows:
            n_no_coords += 1
        clusters.extend(rows)

    if genus is not None:
        LOG.info("genus filter %r: kept %d, skipped %d",
                 genus, len(entries) - n_filtered_genus, n_filtered_genus)
    if n_no_coords:
        LOG.warning("%d entries had no locus with usable coordinates", n_no_coords)

    n = write_ground_truth_tsv(output_path, clusters)
    LOG.info("wrote %d cluster rows → %s", n, output_path)
    return n


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--inspect", type=Path, metavar="DIR",
                   help="print the structure of a few entries and exit "
                        "(use this to verify the schema before building)")
    p.add_argument("--input-dir", type=Path,
                   help="directory of MiBIG per-entry JSON files")
    p.add_argument("--output", type=Path,
                   default=RAW_DIR / "mibig_ground_truth.tsv",
                   help="output TSV (default: %(default)s)")
    p.add_argument("--genus", type=str, default=None,
                   help="keep only entries whose taxonomy name contains this "
                        "string (e.g. 'Streptomyces')")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s  %(levelname)-5s  %(message)s",
        datefmt="%H:%M:%S",
    )

    if args.inspect is not None:
        inspect(args.inspect)
        return

    if args.input_dir is None:
        p.error("either --inspect DIR or --input-dir DIR is required")
    build_ground_truth(args.input_dir, args.output, genus=args.genus)


if __name__ == "__main__":
    main()
