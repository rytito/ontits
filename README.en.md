# ontits

> 🌐 [Español](README.md) (primary version) · **English**

Nextflow pipeline for Oxford Nanopore sequencing of **Illumina-style
dual-indexed fungal ITS amplicons**: from raw POD5 signal to a `phyloseq`
object in a single command.

Derived from [ont16s](../ont_16S_repo/) and built for **Kit 14 chemistry**
(R10.4.1 / FLO-MIN114) with basecalling **disabled** on the sequencer, run on
an HPC cluster.

```
POD5  ──▶  Dorado basecall      (GPU, sup@v5.2.0, --no-trim)
      ──▶  samtools fastq       (whole BAM -> one FASTQ)
      ──▶  cutadapt demux       (dual-index Golay 4F x 4R = 16 samples,
                                 both barcodes required, trims
                                 adapters + barcodes + primers)
      ──▶  chopper              (q≥10, 150–600 bp; ~300 bp amplicon)
      ──▶  seqkit sample        (depth normalisation)
      ──▶  Emu abundance        (UNITE fungi; optional UNITE all)
      ──▶  Emu combine          (taxa × samples tables)
      ──▶  phyloseq             (.rds + QC tables + plots)
```

## The amplicon scheme

The library does **not** use an ONT barcoding kit. It uses the Toju *et al.*
2012 primers in Earth-Microbiome-style dual-index constructs, sequenced on
Nanopore:

```
5'─[P5 + Read1 site (58 nt)]─[Golay F (12 nt)]─[ITS1-F_KYO2 (21 nt)]─
      ...ITS insert...
   ─[rc ITS2_KYO2 (18 nt)]─[rc Golay R (12 nt)]─[rc Read2 site + P7 (58 nt)]─3'
```

- **F primer:** ITS1-F_KYO2 `TAGAGGAAGTAAAAGTCGTAA`
- **R primer:** ITS2_KYO2 `TTYRCTRCGTTCTTCATC`
- **Golay barcodes (12 nt), the same 4 on both ends:**
  `AGCCTTCGTCGC` (1), `TCCATACCGGAA` (2), `AGCCCTGCTACA` (3), `CCTAACGGTCCA` (4)
- **16 combinations** = 16 samples, defined in
  [`assets/barcodes.tsv`](assets/barcodes.tsv). In the original run `F4R3` is
  the mock community and `F4R4` the negative control.

That is why `dorado demux` does not apply here: demultiplexing runs with
**cutadapt** in two rounds. Round 1 orients each read (`--revcomp`) and
assigns it by `barcodeF+primerF`, trimming the 5' end; round 2 assigns it by
`rc(barcodeR+primerR)`, trimming the 3' end. Because the same 4 Golay codes
sit on both ends, a barcode alone cannot identify which end it is — the
adjacent primer can, so barcode and primer are matched as a single adapter.
Only reads with **both** barcodes are accepted (the analogue of
`--barcode-both-ends`); combinations absent from the sheet are discarded as
index hopping.

## Why Emu + UNITE

- **Emu** resolves ambiguous reads via expectation-maximisation over
  full-length alignments — the same reasoning as for 16S applies to ITS,
  where reference databases are full of near-identical sequences between
  sibling species.
- **UNITE** is the standard reference for ITS metabarcoding. Emu ships two
  pre-built versions: `unite-fungi` (fungi only — the primary database) and
  `unite-all` (all eukaryotes — a co-amplification check: KYO2 also amplifies
  plant ITS, and against `unite-fungi` those reads would be forced onto the
  closest fungus).
- **Flat TSVs and a phyloseq object**, not a report you have to scrape.

## Quick start

```bash
# 1. Environments and databases (once per site)
bash bin/setup_environments.sh /shared/conda_envs
micromamba activate /shared/conda_envs/emu
bash bin/setup_databases.sh /shared/databases unite-fungi

# 2. Site configuration
cp conf/site.yml.example conf/site.yml
$EDITOR conf/site.yml          # paths, account, databases

# 3. Run
nextflow run main.nf \
    -profile slurm,conda \
    -params-file conf/site.yml \
    --step all \
    --pod5_dir /project/myrun/pod5 \
    --metadata metadata.tsv \
    --outdir results
```

**Run it inside `tmux`.** Nextflow runs in the foreground and dies with your
SSH session — `-resume` recovers completed work, but the in-flight stage is
lost.

## Stages

Selected with `--step`. Nextflow 26.04 strict syntax removed `-entry`.

| Starting point | Command |
|---|---|
| POD5 (full run) | `--step all --pod5_dir <dir>` |
| An existing `calls.bam` | `--step from_bam --calls_bam calls.bam` |
| Per-sample FASTQ, already trimmed | `--step from_fastq --fastq_dir <dir>` |

## Outputs

```
results/
├── 01_basecall/        tool versions
├── 02_demux/           demux_stats.tsv (reads per sample, unassigned,
│                       invalid combinations) + cutadapt logs
├── 03_filtered/        per-sample FASTQ, barcode/primer-free, filtered
├── 04_subsampled/      depth-normalised FASTQ
├── 05_emu/<db>/        per-sample abundance tables
├── 06_tables/          <db>_emu-combined-{species,genus,family}[-counts].tsv
├── 07_phyloseq/<db>/   phyloseq.rds, alpha diversity, ordination, plots
├── qc/                 read statistics, reads per sample
└── pipeline_info/      timeline, report, trace, DAG
```

Check `02_demux/demux_stats.tsv` first: if a large fraction lands in
`unassigned_*`, lower `--demux_min_overlap` or cautiously raise
`--demux_error_rate` — and reads under `invalid_*` are genuine index hopping.

```r
ps <- readRDS("results/07_phyloseq/unite/phyloseq.rds")
```

## Metadata

A TSV with a `SampleID` column matching the barcode-sheet names
(`F1R1` … `F4R4`). See `assets/metadata_template.tsv`. Every other column
becomes a sample variable; categorical ones are automatically tested with
PERMANOVA.

**Use the mock community (`F4R3`)** to validate the run: it is the
ZymoBIOMICS Microbial Community Standard (Zymo D6300), and with ITS primers
only its two yeasts amplify — the expected result is *Saccharomyces
cerevisiae* + *Cryptococcus neoformans* and nothing else (see
docs/USAGE.en.md). **The blank (`F4R4`)** should fall below `--min_reads`
and be excluded; if it does not, investigate contamination before
interpreting anything.

## Requirements

- Nextflow **≥ 26.04** and Java 11+
- An NVIDIA GPU with tensor cores. **Dorado does not support P100/GP100.**
- Dorado ≥ 2.0 (static binary; deliberately not containerised)
- conda/micromamba, or Apptainer
- cutadapt ≥ 4 (installed by `setup_environments.sh` into the `ont` env)

## Documentation

- [`docs/USAGE.md`](docs/USAGE.md) — parameters, profiles, cluster setup
- [`docs/GOTCHAS.md`](docs/GOTCHAS.md) — failure modes we hit, and how to
  detect them

## How to cite

Cite the underlying tools: Dorado (Oxford Nanopore), cutadapt (Martin 2011),
Emu (Curry *et al.* 2022, *Nature Methods*), minimap2 (Li 2018), chopper,
seqkit, phyloseq (McMurdie & Holmes 2013), the ITS1-F_KYO2 / ITS2_KYO2 primers
(Toju *et al.* 2012, *PLoS ONE*), plus **UNITE with its version number**
(Nilsson *et al.* 2019; Emu's pre-built databases use approximately the 2021
release — verify and report it in your methods).

## Funding and credits

Developed within the **MicroAndes** project — *Harnessing microbes to
strengthen Andean communities of fermented-food micro-producers* — a
**VLIR-UOS Short Initiative** (KU Leuven project
[3M250529](https://research.kuleuven.be/portal/nl/project/3M250529), funding
reference PE2025SIN468A101, 2025–2027). Coordinated by KU Leuven (Molecular
Bacteriology Lab, Rega Institute) with Universidad
Nacional del Centro del Perú (UNCP) and Andean micro-producer communities.

## License

MIT — see [`LICENSE`](LICENSE).
