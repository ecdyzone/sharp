#!/usr/bin/env python3
"""Index a directory of NCBI assembly dumps into the files the pipeline needs.

The input is a genome database laid out the way `datasets`/the NCBI FTP mirror
leaves it: one directory per assembly, each holding `<dir>_genomic.fna` (plus
`.gbff`/`.gff`/`protein.faa`, which this script does not read).

    <db>/GCF_000009765.2_ASM976v2/GCF_000009765.2_ASM976v2_genomic.fna
    <db>/GCF_041549525.1_P_aurescens_B2879/..._genomic.fna

Only the FASTA is parsed, and only for its headers and sequence lengths. That
is deliberate: the `.gbff` files are 2-5x larger and carry nothing this needs.

Why this exists: the rest of the benchmark pipeline is keyed by **contig
accession**, while a genome database is keyed by **assembly**. One assembly
holds many contigs, so neither key alone can join a prediction back to its
genome. The FASTA headers carry the versioned nucleotide accession verbatim
(`>NC_003155.5 Streptomyces avermitilis ...`), which is exactly the name every
tool will echo back in its output — so this script reads the bridge straight
out of the data rather than needing an external index.

Four outputs, each with one job:

    genomes.tsv          assembly, contig, length, description   (per contig)
    assemblies.tsv       assembly, n_contigs, total_bp, organism, fasta
    assemblies.txt       one assembly per line  -> the job-array line list
    analyzed_contigs.txt one contig per line    -> the --contigs scope file

`--link-dir` additionally builds a symlink farm, `<link-dir>/<assembly>.fasta`
-> the real `.fna`. That is what lets the existing array scripts run over this
database unchanged: they expect `<GENOME_DIR>/<key>.fasta`, and the symlink
also makes each tool name its output after the assembly, so the output pool
stays keyed the same way as the manifest.

Usage:
    # Verify the layout is understood before indexing 2,500 genomes
    python scripts/build_genome_manifest.py --db-dir /path/to/db --inspect

    # Index, and stage the symlinks the array scripts read
    python scripts/build_genome_manifest.py \\
        --db-dir /path/to/actinomycetota_db \\
        --output-dir data/interim/actino_db \\
        --link-dir data/raw/actino_db
"""
from __future__ import annotations

import argparse
import gzip
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from sharp.config import INTERIM_DIR  # noqa: E402

LOG = logging.getLogger("build_genome_manifest")

# ═══════════════════════════ layout assumptions ════════════════════════════
# Everything this script assumes about the database layout lives here. If a
# differently-shaped dump has to be indexed, this is the block to change.

# Preference order: the NCBI `_genomic.fna` naming first, then looser fallbacks
# so a hand-assembled directory still works.
FASTA_PATTERNS = ("*_genomic.fna", "*_genomic.fna.gz", "*.fna", "*.fna.gz",
                  "*.fasta", "*.fa")

# `GCF_000009765.2_ASM976v2` -> `GCF_000009765.2`. The assembly *accession* is
# the stable, universally-quoted key; the `_ASM976v2` suffix is the assembly
# *name* and only exists in the directory listing. Keying on the accession is
# what lets a collaborator running a third tool join results without having to
# reproduce our directory names.
ACCESSION_RE = re.compile(r"^(GC[AF]_\d+\.\d+)")

# Trailing noise in a FASTA description, stripped to get something usable as an
# organism label. Heuristic and only used for grouping/reporting — nothing
# downstream depends on it being exactly right.
ORGANISM_NOISE_RE = re.compile(
    r"\s*,?\s*(complete (genome|sequence|cds)|chromosome|plasmid\b.*|"
    r"scaffold\b.*|contig\b.*|whole genome shotgun sequence)\s*$",
    re.IGNORECASE,
)
# ═══════════════════════════════════════════════════════════════════════════


@dataclass(frozen=True)
class ContigRecord:
    """One FASTA record, tied back to the assembly that contains it."""
    assembly: str
    contig: str
    length: int
    description: str


@dataclass(frozen=True)
class AssemblyRecord:
    """One assembly, summarised over its contigs."""
    assembly: str
    directory: str
    fasta: Path
    n_contigs: int
    total_bp: int
    organism: str


def assembly_key(directory_name: str) -> str:
    """`GCF_000009765.2_ASM976v2` -> `GCF_000009765.2`.

    Falls back to the whole directory name for a dump that is not named after
    an NCBI accession, so a non-standard directory is indexed rather than
    silently dropped.
    """
    m = ACCESSION_RE.match(directory_name)
    return m.group(1) if m else directory_name


def find_fasta(directory: Path) -> Path | None:
    """The genomic FASTA inside one assembly directory, or None."""
    for pattern in FASTA_PATTERNS:
        hits = sorted(directory.glob(pattern))
        # `_protein.faa` is excluded by the patterns; a `_cds_from_genomic`
        # companion is not, so prefer the shortest name, which is the plain
        # genomic FASTA.
        hits = [h for h in hits if "cds_from" not in h.name
                and "rna_from" not in h.name]
        if hits:
            return min(hits, key=lambda p: len(p.name))
    return None


def _open_text(path: Path):
    """Open a FASTA whether or not it is gzipped."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_fasta_index(path: Path, assembly: str) -> list[ContigRecord]:
    """Header + sequence length for every record, without holding sequence.

    Lengths are counted as they stream past, so a 12 Mb genome costs one line
    of memory rather than 12 MB.
    """
    records: list[ContigRecord] = []
    contig: str | None = None
    description = ""
    length = 0

    def flush() -> None:
        if contig is not None:
            records.append(ContigRecord(assembly, contig, length, description))

    with _open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                header = line[1:].strip()
                parts = header.split(None, 1)
                contig = parts[0] if parts else ""
                description = parts[1] if len(parts) > 1 else ""
                length = 0
            else:
                length += len(line.strip())
    flush()
    return records


def organism_of(records: list[ContigRecord]) -> str:
    """A readable organism label, taken from the longest contig's description.

    The longest record is the chromosome in an assembly that has one, which
    carries a better label than a plasmid or a short scaffold.
    """
    if not records:
        return ""
    longest = max(records, key=lambda r: r.length)
    label = longest.description
    # Strip trailing descriptors repeatedly: "…, chromosome, complete genome".
    for _ in range(3):
        stripped = ORGANISM_NOISE_RE.sub("", label).rstrip(" ,")
        if stripped == label:
            break
        label = stripped
    return label


def assembly_dirs(db_dir: Path) -> list[Path]:
    """Every subdirectory of the database root, sorted."""
    if not db_dir.is_dir():
        raise SystemExit(f"not a directory: {db_dir}")
    return sorted(d for d in db_dir.iterdir() if d.is_dir())


def index_database(
    db_dir: Path, limit: int | None = None
) -> tuple[list[AssemblyRecord], list[ContigRecord], list[str]]:
    """Walk the database, returning (assemblies, contigs, problems)."""
    assemblies: list[AssemblyRecord] = []
    contigs: list[ContigRecord] = []
    problems: list[str] = []
    seen: dict[str, str] = {}

    dirs = assembly_dirs(db_dir)
    if limit is not None:
        dirs = dirs[:limit]

    for i, d in enumerate(dirs, 1):
        key = assembly_key(d.name)
        if key in seen:
            problems.append(
                f"{d.name}: duplicate assembly key {key} "
                f"(already seen as {seen[key]}) — skipped"
            )
            continue

        fasta = find_fasta(d)
        if fasta is None:
            problems.append(f"{d.name}: no genomic FASTA found — skipped")
            continue

        recs = parse_fasta_index(fasta, key)
        if not recs:
            problems.append(f"{d.name}: {fasta.name} holds no records — skipped")
            continue

        seen[key] = d.name
        contigs.extend(recs)
        assemblies.append(AssemblyRecord(
            assembly=key,
            directory=d.name,
            fasta=fasta.resolve(),
            n_contigs=len(recs),
            total_bp=sum(r.length for r in recs),
            organism=organism_of(recs),
        ))
        if i % 200 == 0:
            LOG.info("indexed %d/%d directories", i, len(dirs))

    return assemblies, contigs, problems


# ══════════════════════════════ writing ════════════════════════════════════

def write_genomes_tsv(path: Path, contigs: list[ContigRecord]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("assembly\tcontig\tlength\tdescription\n")
        for c in contigs:
            fh.write(f"{c.assembly}\t{c.contig}\t{c.length}\t{c.description}\n")
    return len(contigs)


def write_assemblies_tsv(path: Path, assemblies: list[AssemblyRecord]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        fh.write("assembly\tn_contigs\ttotal_bp\torganism\tfasta\n")
        for a in assemblies:
            fh.write(f"{a.assembly}\t{a.n_contigs}\t{a.total_bp}\t"
                     f"{a.organism}\t{a.fasta}\n")
    return len(assemblies)


def write_lines(path: Path, values: list[str]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("".join(f"{v}\n" for v in values))
    return len(values)


def load_genomes_tsv(path: Path) -> dict[str, list[str]]:
    """Read genomes.tsv back as {assembly: [contig, ...]}.

    Shared with merge_predictions.py and summarize_predictions.py, so the
    manifest format is parsed in exactly one place.
    """
    out: dict[str, list[str]] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_asm, i_ctg = header.index("assembly"), header.index("contig")
        except ValueError as e:
            raise ValueError(
                f"{path} is not a genomes.tsv manifest "
                f"(need 'assembly' and 'contig' columns, got {header})"
            ) from e
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out.setdefault(f[i_asm], []).append(f[i_ctg])
    return out


def load_contig_index(path: Path) -> dict[str, tuple[str, int]]:
    """Read genomes.tsv back as {contig: (assembly, length)}.

    This is the join every downstream comparison needs: a tool reports a
    contig, and the contig has to become an assembly (to group by genome) and a
    length (to normalise counts per Mb).
    """
    out: dict[str, tuple[str, int]] = {}
    with path.open() as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_asm = header.index("assembly")
            i_ctg = header.index("contig")
            i_len = header.index("length")
        except ValueError as e:
            raise ValueError(
                f"{path} is not a genomes.tsv manifest (need 'assembly', "
                f"'contig' and 'length' columns, got {header})"
            ) from e
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            out[f[i_ctg]] = (f[i_asm], int(f[i_len]))
    return out


def link_genomes(link_dir: Path, assemblies: list[AssemblyRecord]) -> tuple[int, list[str]]:
    """Stage `<link_dir>/<assembly>.fasta` symlinks into the database.

    Symlinks rather than copies: the database is tens of gigabytes and already
    lives on the machine that will run the tools.
    """
    link_dir.mkdir(parents=True, exist_ok=True)
    skipped: list[str] = []
    n = 0
    for a in assemblies:
        if a.fasta.suffix == ".gz":
            # Not all tools read gzipped FASTA, and silently linking one would
            # fail 2,000 array tasks in. Report instead.
            skipped.append(f"{a.assembly}: {a.fasta.name} is gzipped — "
                           f"gunzip it before running the tools")
            continue
        link = link_dir / f"{a.assembly}.fasta"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(a.fasta.resolve())
        n += 1
    return n, skipped


# ══════════════════════════════ inspection ═════════════════════════════════

def inspect(db_dir: Path, limit: int) -> None:
    """Print what the layout looks like, without writing anything."""
    dirs = assembly_dirs(db_dir)
    print(f"\n{'='*72}\nDATABASE LAYOUT — {db_dir}\n{'='*72}")
    print(f"{len(dirs)} assembly directories")
    if not dirs:
        return

    print(f"\nFiles in the first directory ({dirs[0].name}):")
    for f in sorted(dirs[0].iterdir()):
        print(f"  {f.name:<60} {f.stat().st_size/1e6:>8.1f} MB")

    assemblies, contigs, problems = index_database(db_dir, limit=limit)
    print(f"\n{'-'*72}\nINDEXED (first {limit})\n{'-'*72}")
    print(f"{'assembly':<20} {'contigs':>7} {'total bp':>12}  organism")
    for a in assemblies:
        print(f"{a.assembly:<20} {a.n_contigs:>7} {a.total_bp:>12,}  {a.organism}")

    print(f"\n{'-'*72}\nCONTIG ACCESSIONS (the join key)\n{'-'*72}")
    for c in contigs[:12]:
        print(f"  {c.assembly:<20} {c.contig:<18} {c.length:>12,}  {c.description[:40]}")
    if len(contigs) > 12:
        print(f"  ... and {len(contigs)-12} more")

    if problems:
        print(f"\n{'-'*72}\nPROBLEMS\n{'-'*72}")
        for p in problems:
            print(f"  {p}")
    print()


# ══════════════════════════════ orchestration ══════════════════════════════

def run(
    db_dir: Path,
    output_dir: Path | None,
    link_dir: Path | None,
    limit: int | None,
    do_inspect: bool,
) -> None:
    if do_inspect:
        inspect(db_dir, limit or 5)
        return

    if output_dir is None:
        raise SystemExit("--output-dir is required unless --inspect is given")

    assemblies, contigs, problems = index_database(db_dir, limit=limit)
    if not assemblies:
        raise SystemExit(f"no usable assembly directories under {db_dir}")

    n_c = write_genomes_tsv(output_dir / "genomes.tsv", contigs)
    n_a = write_assemblies_tsv(output_dir / "assemblies.tsv", assemblies)
    write_lines(output_dir / "assemblies.txt", [a.assembly for a in assemblies])
    write_lines(output_dir / "analyzed_contigs.txt", [c.contig for c in contigs])

    total_bp = sum(a.total_bp for a in assemblies)
    LOG.info("indexed %d assemblies / %d contigs / %.1f Gb", n_a, n_c, total_bp / 1e9)
    LOG.info("wrote %s/{genomes.tsv,assemblies.tsv,assemblies.txt,"
             "analyzed_contigs.txt}", output_dir)

    # A fragmented assembly cannot yield an intact BGC for any tool, so the
    # spread is worth knowing before the numbers are interpreted.
    fragmented = [a for a in assemblies if a.n_contigs > 50]
    if fragmented:
        LOG.warning("%d/%d assemblies have >50 contigs (most fragmented: %s, %d) "
                    "— draft assemblies split clusters across contig boundaries "
                    "and depress every tool's counts",
                    len(fragmented), len(assemblies),
                    max(fragmented, key=lambda a: a.n_contigs).assembly,
                    max(a.n_contigs for a in fragmented))

    if problems:
        LOG.warning("%d directory/-ies were skipped:", len(problems))
        for p in problems[:20]:
            LOG.warning("  %s", p)
        if len(problems) > 20:
            LOG.warning("  ... and %d more", len(problems) - 20)

    if link_dir is not None:
        n_links, skipped = link_genomes(link_dir, assemblies)
        LOG.info("staged %d symlinks → %s/<assembly>.fasta", n_links, link_dir)
        for s in skipped[:20]:
            LOG.warning("  %s", s)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--db-dir", type=Path, required=True,
                   help="database root holding one directory per assembly")
    p.add_argument("--output-dir", type=Path, default=None,
                   help="where the four manifest files go "
                        f"(e.g. {INTERIM_DIR / 'actino_db'})")
    p.add_argument("--link-dir", type=Path, default=None,
                   help="also stage <link-dir>/<assembly>.fasta symlinks for "
                        "the array scripts to read")
    p.add_argument("--limit", type=int, default=None,
                   help="index only the first N directories (dev/pilot)")
    p.add_argument("--inspect", action="store_true",
                   help="print the database layout and write nothing")
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
        db_dir=args.db_dir,
        output_dir=args.output_dir,
        link_dir=args.link_dir,
        limit=args.limit,
        do_inspect=args.inspect,
    )


if __name__ == "__main__":
    main()
