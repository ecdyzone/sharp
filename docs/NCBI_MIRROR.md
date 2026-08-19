# The local NCBI mirror on davinci

Findings from inspecting the university cluster's NCBI FTP mirror (2026-08-19),
recorded so the question "could we read genomes from the mirror instead of
downloading them?" does not have to be re-investigated from scratch.

**Conclusion: not for this benchmark.** The mirror is complete and well-formed;
it is keyed at the wrong level. See [Why we download anyway](#why-we-download-anyway).
Revisit if the benchmark scales past ~1000 genomes.

---

## What is actually there

`$DATABASES` on davinci is **`/scratch/global/databases`** — not `/databases`,
which is what a snippet circulating internally hardcodes. Anything reading this
mirror must use the environment variable, not a literal path.

```
$DATABASES/genomes/
├── all/
│   ├── GCA/                 ← GenBank assemblies    (present)
│   └── GCF/                 ← RefSeq assemblies     (present)
├── ASSEMBLY_REPORTS/
│   ├── assembly_summary_genbank.txt      (present)
│   ├── assembly_summary_refseq.txt       (present)
│   ├── assembly_summary_*_historical.txt (suppressed/replaced assemblies)
│   ├── ftp_mod_times.txt.gz              (per-assembly mirror timestamps)
│   └── ... ANI / CheckM / FCS QC reports
└── TARGET/                  ← rRNA reference sets only (5S/18S/23S/28S/ITS),
                               NOT a sequence-accession index
```

Verified facts:

| Question | Answer |
|---|---|
| `$DATABASES` | `/scratch/global/databases` |
| Tree layout | `genomes/all/{GCA,GCF}/NNN/NNN/NNN/<asm>/` — matches the FTP URL rewrite |
| Per-assembly `*_assembly_report.txt` | **Present**, on both the GCA and GCF sides |
| `assembly_summary_genbank.txt` | **Present** — so no GenBank→RefSeq hop is needed at the *assembly* level |
| Mirror freshness | Both summaries stamped **2025-10-26** (~10 months stale as of 2026-08-19) |
| Scale | 12,852 *Streptomyces* rows in `assembly_summary_genbank.txt` |

Both `GCA_000203835.1` and `GCF_000203835.1` (*S. coelicolor*) carry the full
file set: `_genomic.fna.gz`, `_genomic.gbff.gz`, `_assembly_report.txt`,
checksums. Nothing is missing from the mirror itself.

## The problem: assembly-keyed vs nucleotide-keyed

GCA/GCF accessions name a **whole assembly**. Our ground truth names
**individual nucleotide records** (`AL645882.2`, `CP002047.1`, …). These are
different levels of the same hierarchy, and the mirror is indexed only by the
upper one.

**Neither `assembly_summary_*.txt` contains nucleotide accessions at all.** They
map assembly → organism → FTP path. So `AL645882.2` cannot be looked up in
either summary table, and no other file in `ASSEMBLY_REPORTS/` provides that
index — a search for `nucl|accession|seq|contig` matches only the RefSeq
assembly summaries and FCS QC reports, and `TARGET/` turned out to be rRNA
reference sets.

The nucleotide→assembly mapping exists **only inside each per-assembly
`*_assembly_report.txt`**, which lists that assembly's sequences with both their
GenBank and RefSeq accessions. Building the reverse index therefore means
reading ~2M of those files across the tree (or, for a *Streptomyces*-only
restriction, ~12.8k of them). That is feasible but it is a real preprocessing
step with its own artifact to keep fresh — not a lookup.

## Two further problems specific to this benchmark

**1. WGS accessions do not name a whole assembly.** 8 of our 50 accessions are
scaffold-level (`KZ*`, `CM*`, `GG*`, `VAWE*`, `RCOL*`, `PTJS*`, `JJOB*`, `DS*`).
The mirror stores the assembly's *complete* FASTA — hundreds of contigs, only
one of which the ground truth names. Using it means subsetting back down to the
named record, which is a second silent-failure mode stacked on the first: get it
wrong and the contig names stop matching the ground truth, and the affected
genomes score zero with nothing in the metrics to indicate why.

**2. Reproducibility.** The mirror is an rsync of a moving target, currently
~10 months stale. A benchmark number must be reproducible from its command line
(see CLAUDE.md, "`.env` is for machine identity, not pipeline behaviour"). A
genome resolved through "whatever the mirror held when it last synced" does not
satisfy that; `accession.version` does.

## Why we download anyway

`scripts/download_benchmark_genomes.sh` fetches by nucleotide accession, so the
FASTA header **is** the name the ground truth uses — and the script verifies
exactly that, per genome. No translation layer, no subsetting, no index to
maintain, and the benchmark pins to `accession.version`.

The cost avoided is ~420 Mb, once, for 50 genomes.

## If we do use the mirror later

The threshold is scale: at thousands of genomes, downloading stops being
reasonable and the preprocessing cost amortizes. It would be **one localized
change** — everything downstream of `data/raw/genomes/<ACC>.fasta`
(arrays → `merge_predictions.py` → `sharp.evaluate`) is indifferent to how the
FASTA arrived. Sketch:

1. Build a nucleotide→assembly index once by walking
   `$DATABASES/genomes/all/GC[AF]/*/*/*/*/*_assembly_report.txt`, restricted to
   *Streptomyces* rows of `assembly_summary_genbank.txt` to keep it to ~12.8k
   files. Cache it as a TSV under `data/interim/`, the same way
   `record_lengths.tsv` caches NCBI lengths.
2. Add a mirror branch to `download_benchmark_genomes.sh`: resolve the
   accession through the index, `zcat` the assembly's `_genomic.fna.gz`, and
   **subset to the named record** for WGS accessions.
3. Keep the existing header verification unconditional — it is what would catch
   a bad translation, and it matters more on the mirror path than on the
   download path.
4. Fall back to `fetch_nuccore` on any index miss, so mirror staleness degrades
   to a download rather than to a missing genome.
