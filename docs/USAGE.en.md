# Usage

> 🌐 [Español](USAGE.md) (primary version) · **English**

Initial site setup (Java, Nextflow, Dorado, shell environment, storage
planning) is identical to ont16s — see the
[ont16s USAGE](../../ont_16S_repo/ont16s/docs/USAGE.en.md) for full detail.
Below is what is specific to ontits.

## Environments and databases

```bash
bash bin/setup_environments.sh /shared/conda_envs      # adds cutadapt to the ont env
micromamba activate /shared/conda_envs/emu
bash bin/setup_databases.sh /shared/databases unite-fungi   # or: unite-all | all
```

Emu's pre-built ITS databases come from UNITE:

| Database | Use it for | Do not use it for |
|---|---|---|
| **unite-fungi** | Primary taxonomy (Fungi only) | Detecting plant co-amplification |
| **unite-all** | Quantifying co-amplification (all eukaryotes) | — |

Running both is the ITS version of ont16s's GTDB/SILVA check: KYO2
co-amplifies plant and other eukaryote ITS; against `unite-fungi` those reads
are forced onto the closest fungus or lost. Both databases are small
(≤100 MB), so running both is cheap.

The pre-builts use the ~2021 UNITE release. For the current version, build
with `emu build-database` from the UNITE general FASTA release.

## Parameters

### Input

| Parameter | Default | Notes |
|---|---|---|
| `--step` | `all` | `all` \| `from_bam` \| `from_fastq` |
| `--pod5_dir` | — | With basecalling disabled, MinKNOW writes `pod5/`, not `pod5_pass/` |
| `--calls_bam` | — | For `--step from_bam` |
| `--fastq_dir` | — | For `--step from_fastq`; files `<SampleID>.fastq`, **already barcode/primer-free** |
| `--metadata` | — | TSV with `SampleID` (`F1R1`…`F4R4`); enables phyloseq |
| `--outdir` | `results` | |

### Basecalling

Identical to ont16s: `--model sup` (→ `sup@v5.2.0` for R10.4.1),
`--max_reads` for tests. `--no-trim` is hard-wired: Dorado's automatic
adapter trimming would eat the barcodes cutadapt needs.

### Demultiplexing (cutadapt dual-index)

| Parameter | Default | Notes |
|---|---|---|
| `--barcodes` | `assets/barcodes.tsv` | Sheet of 16 combinations (SampleID, fwd_id, rev_id, fwd_barcode, rev_barcode) |
| `--fwd_primer` | `TAGAGGAAGTAAAAGTCGTAA` | ITS1-F_KYO2 |
| `--rev_primer` | `TTYRCTRCGTTCTTCATC` | ITS2_KYO2 (IUPAC: Y, R) |
| `--demux_error_rate` | `0.15` | ~4 errors tolerated over the 30–33 nt barcode+primer adapter |
| `--demux_min_overlap` | `25` | Requires seeing the barcode **and** primer, not a fragment |

Two rounds: (1) `--revcomp` orients the read and assigns the F barcode
(adjacent to its primer, matched as a single adapter — with the same 4 Golay
codes on both ends, a barcode alone cannot identify the end); (2) assigns the
R barcode at the 3' end. Only reads with **both** barcodes are accepted;
combinations absent from the sheet are discarded as index hopping and counted
as `invalid_*` in `demux_stats.tsv`.

If `unassigned_*` exceeds ~25%: check run quality first; then lower
`--demux_min_overlap` (e.g. to 20) before raising the error rate — the Golay
codes differ at ≥5 positions, but that margin erodes quickly.

### Filtering and depth

| Parameter | Default | Notes |
|---|---|---|
| `--min_qscore` | `10` | |
| `--min_len` / `--max_len` | `150` / `600` | The KYO2 amplicon is ~300 bp (insert ~260 bp once primers/barcodes are trimmed); the window leaves room for genuine ITS1 length variation — do not narrow it without checking `qc/read_stats.tsv` |
| `--subsample` | `20000` | `0` disables |
| `--seed` | `42` | Reproducible subsampling |
| `--min_reads` | `1000` | Lower than 16S; the negative control (`F4R4`) should fall below it |

**Unlike 16S**, ITS1 length varies biologically across taxa. A narrow window
biases composition against long- or short-ITS taxa — hence the deliberately
wide range, which only removes fragments and concatemers.

### Taxonomy

| Parameter | Notes |
|---|---|
| `--emu_dbs` | `name:path[,...]` — e.g. `'unite:/db/emu_unite_fungi'` |

## Profiles

Identical to ont16s: `standard`, `slurm`, `vsc_wice`, `vsc_genius`, `conda`,
`apptainer`, `test`, `debug`. Combine with commas: `-profile vsc_wice,conda`.

## Run quality control

1. **`results/02_demux/demux_stats.tsv`** — first. Reads per sample,
   `unassigned_*` (missing one of the two barcodes) and `invalid_*` (unused
   combination = measured index hopping).
2. **Mock community (`F4R3`)** — expected vs observed composition validates
   primers, demux and database at once.
3. **Negative control (`F4R4`)** — should fall below `--min_reads` and be
   excluded from the phyloseq object. If not: contamination.
4. **`07_phyloseq/<db>/uncharacterised_fraction.tsv`** — per-sample
   UNITE "unidentified"/SH fraction. In understudied systems this fraction is
   a result, not an error.

## Safe execution and resume

Same as ont16s: run inside `tmux`, never Ctrl+Z, resume with `-resume`.
After a successful basecall and a later failure:

```bash
nextflow run main.nf ... --step from_bam --calls_bam results/01_basecall/calls.bam -resume
```
