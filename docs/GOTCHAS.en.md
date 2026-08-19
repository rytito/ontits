# Gotchas

> **Note (ontits):** inherited from ont16s. Everything about disk quotas, Dorado, SLURM/VSC, Nextflow, Emu, tmux and data transfer applies as-is. Sections about `dorado demux` / `--kit-name`, GTDB/SILVA and 16S species resolution do NOT apply here: ontits demultiplexes with cutadapt (see README) and profiles against UNITE. The lessons on unused barcode combinations as a cross-talk measure and on low-depth samples apply doubly to ITS.

> 🌐 [Español (versión principal)](GOTCHAS.md) · **English**

Failure modes we have actually hit, with how to spot each one. Several report
success while silently losing data.

---

## Silent data loss

### Disk quota truncates BAMs, and Dorado still exits 0

**Symptom:** `EOF marker is absent. The input is probably truncated`, while
`sacct` shows `COMPLETED 0:0`.

**Cause:** VSC scratch is **100 GB**, not 500. With ~86 GB of POD5 plus
`calls.bam` plus demux output, writes hit the wall. Dorado does not propagate
the write failure as a non-zero exit status.

**Detect:**
```bash
myquota
samtools quickcheck -u -v calls.bam
find demux -name '*.bam' | xargs samtools quickcheck -u -v && echo "ALL OK"
```

The `-u` matters: without it, `quickcheck` reports `had no targets in header`,
which is **normal** for unaligned BAM. Only "missing EOF block" means
truncation.

**Fix:** keep POD5 on project storage and point `--pod5_dir` there. The
pipeline sets `stageInMode = 'symlink'` so Nextflow does not duplicate it.

### Reading output while a job is still writing

The same truncation errors appear if you count reads in `demux/` while
`dorado demux` is still running. Check `squeue -M wice -u $USER` before
concluding anything is broken.

---

## Dorado

### `--kit-name` on the basecaller does not split reads

It only *tags* them. `dorado demux` is still required to produce per-sample
files. The most common misconception about the tool.

### `--no-trim` is mandatory when demuxing later

Without it the basecaller strips adapters and barcodes, and demux finds
nothing — everything lands in `unclassified`.

### Output is a nested MinKNOW tree

```
demux/<experiment>/<sample_id>/<run_id>/bam_pass/barcodeNN/*.bam
```

Rebuilt from POD5 metadata. When MinKNOW had no sample ID set, the middle
level is literally `unknown`. **Never glob a fixed path** — traverse with
`find`. Any script using `demux/*.fastq` silently matches nothing.

### P100 GPUs are unsupported

Dorado's requirements exclude P100/GP100. On Genius this is often the emptiest
partition and therefore the tempting one. Use `gpu_v100`, or wICE.

### `ReqNodeNotAvail, Reserved for maintenance` with no reservation

`scontrol show reservation` may report nothing while the job still refuses to
schedule — the nodes are drained.

```bash
sinfo -M wice -p gpu_h100 -o "%P %a %D %t %N"
sinfo -M wice -R
```

`maint` and a trailing `$` mean out of service; `mix`/`alloc` just mean busy.
`SCHEDNODES (null)` with `START_TIME N/A` means no node could ever satisfy the
request.

---

## SLURM on VSC

### `--nodes` must be explicit on GPU submissions

```
sbatch: error: [ERROR] Please specify the number of nodes (using --nodes)
```

Nextflow derives `--nodes` from `cpus` for ordinary processes, but a
`clusterOptions` override replaces the whole string. The GPU label must carry
`--nodes=1 --ntasks=1 --gpus-per-node=1` explicitly.

### Per-GPU core ceilings

Exceeding these gets the job rejected outright.

| Cluster | Partition | Max cores/GPU | Max mem/GPU |
|---|---|---|---|
| wICE | `gpu_h100` | 16 | 187200 MiB |
| wICE | `gpu_a100` | 18 | 126000 MiB |
| Genius | `gpu_v100` | 4 | 84000 MiB |

### `-M` / `--clusters` is mandatory

Plain `squeue` defaults to the login node's cluster and looks empty even when
jobs are running. Always `squeue -M wice -u $USER`.

### `--parsable` returns `jobid;cluster`

So `--dependency=afterok:$JID` fails with "Job dependency problem":

```bash
RAW=$(sbatch --parsable job.slurm)
JID=${RAW%%;*}
```

### `#!/bin/bash -l` is required

Without `-l` the login shell never initialises, the cluster module does not
load, and `PATH` additions from `~/.bashrc` are absent.

### Backgrounded `srun` does not survive logout

An allocation suspended with Ctrl+Z and resumed with `bg` dies when the SSH
session drops. Use `sbatch`.

### Pasting a block while `srun` queues runs it on the login node

The `srun` waits; every following line buffers in the login shell and executes
there. This is how Dorado ends up reporting "no CUDA devices available" on a
login node. Submit one line, wait for the prompt to change, confirm with
`hostname` and `nvidia-smi`, then continue.

---

## Nextflow

### `-entry` removed in 26.04

```
ERROR ~ The `-entry` option is not supported with the strict parser
```

The strict syntax (parser v2) became the default. This pipeline uses `--step`
instead — a pipeline param, double dash:

```bash
nextflow run main.nf --step from_fastq --fastq_dir in
```

### Other strict-syntax rejections

All four hit during migration:

| Construct | Replacement |
|---|---|
| `exit 1, "msg"` | `error "msg"` |
| `switch`/`case` | `if` / `else if` chain |
| `workflow.onComplete` | delete it; Nextflow prints its own summary |
| top-level statements | move inside a workflow or function |

Bash `exit 1` **inside** a `script:` block is untouched by this — only Groovy
code is affected.

### `publishDir` with input variables needs a closure

```groovy
publishDir "${params.outdir}/${db_name}"      // No such variable: db_name
publishDir { "${params.outdir}/${db_name}" }  // correct
```

Directives are evaluated before input variables bind. Same applies to
`saveAs`, where the closure parameter should be named explicitly
(`{ fn -> ... }`) rather than relying on implicit `it`.

### Command-line params arrive as Strings

```
Cannot compare java.lang.String with value '5000' and java.lang.Integer
```

`--subsample 5000` is a String; `params.subsample > 0` then fails. Use
`params.subsample.toString().toInteger()`. Adding `nf-schema` with typed
declarations would handle this properly.

### Never set `PATH` in an `env {}` block

Nextflow emits `env` exports **after** `beforeScript`, so a `PATH` assignment
there silently overwrites anything `beforeScript` set. Symptom: `.command.run`
contains two `export PATH=` lines, the second lacking your tools. Put tool
paths in `process.beforeScript` only.

### conda and micromamba are shell *functions* on VSC

```bash
which micromamba   # prints a function body, not a path
```

Nextflow processes never source `~/.bashrc`, so the `conda` directive cannot
find them. Put each environment's `bin/` on `PATH` via `beforeScript` instead.

### `withName: '!EMU_.*'` is not valid

Negated name selectors silently match nothing. Set a process-level default and
override for the specific case — later, more specific selectors win.

### Do not Ctrl+Z a running pipeline

It leaves the session lock held:

```
ERROR ~ Unable to acquire lock on session with ID ...
```

Fix with `jobs` then `kill %1`. Ctrl+C once instead — Nextflow cancels
cleanly and releases the lock.

---

## Emu

### SILVA needs ~48 GB of RAM

SILVA 138.2 has **2.2 million reference sequences (3.2 Gbp)**. The minimap2
index alone exceeds 16 GB, and a 16 GB allocation gets OOM-killed:

```
Command 'minimap2 -ax map-ont ...' died with <Signals.SIGKILL: 9>
```

The `cpu_med` label requests 48 GB and doubles on retry. GTDB (77k sequences)
runs comfortably in 16 GB. Including SILVA roughly triples the cost of a run —
worth it for the chloroplast check, unnecessary for routine profiling.

### Fractional counts

EM outputs estimates like `39879.125`. Round before DESeq2 or
`rarefy_even_depth()`, both of which need integers.

### `--keep-counts` is required for `combine-outputs --counts`

Without it you get relative abundances only.

### Runtime scales steeply with depth

90k reads against SILVA took >20 min/sample and 284 EM iterations, versus 32
against the default database. Subsample to 20k — relative abundances are
stable well below that.

### Python version

Emu's README targets Python 3.7; bioconda installs it under 3.13. It works,
but export `PYTHONNOUSERSITE=1` so `~/.local` packages do not shadow the
environment. If imports fail, rebuild pinned to 3.10.

---

## Databases

### GTDB species are binomials — do not split headers on whitespace

```
>RS_GCF_000566285.1 d__Bacteria;...;s__Flavobacterium fluviatile [ssu_len=1510]
```

A naive `header.split()` truncates the species to `Flavobacterium` and
collapses every species in a genus into one taxon. **Symptom:** far fewer
unique taxa than sequences — we saw 20,910 from 75,602. Extract from `d__` up
to the first ` [`. `bin/gtdb_to_emu.py` does this correctly.

**Verify after building:**
```bash
head -2 emu_input/taxonomy.tsv | awk -F'\t' '{print $2}'
# expect "Flavobacterium fluviatile", not "Flavobacterium"
```

### `column -t` misrenders the taxonomy table

It splits on any whitespace, so the space inside a binomial looks like a column
break and every rank appears shifted. The file is fine:

```bash
head -3 taxonomy.tsv | awk -F'\t' '{for(i=1;i<=NF;i++) printf "[%s] ", $i; print ""}'
```

### Use `ssu_reps`, not `ssu_all`

`ssu_all` holds every SSU copy from 732k genomes. `ssu_reps` is one per species
cluster — ~77k sequences, comparable to Emu's default.

### GTDB genus suffixes are real

`Clostridium_S`, `Clostridium_B`, `Clostridium_AM` are **different genera**,
appended where the NCBI genus is polyphyletic. Never strip them with a regex —
it silently merges distinct lineages. Flag them in your methods so reviewers do
not read them as typos.

### GTDB is prokaryote-only

No chloroplast or mitochondrial references. Plant reads map to the nearest
cyanobacterium and get a plausible-looking name. Use SILVA to quantify that
fraction.

---

## Interpretation

### Low-depth samples show *maximum* apparent diversity

A sample with 669 reads gave 27 taxa above 1% — the highest in the run — while
a healthy sample gave 7. Maximum diversity at minimum depth is the signature of
noise. Drop samples below ~5000 reads.

### Unused barcodes are your only cross-talk measurement

Filling all 24 slots leaves no negative control and no detection floor. Leave
one empty, or include an extraction blank.

### Closely related species cannot be separated by 16S

The *L. casei* group and the solventogenic clostridia have near-identical 16S
genes. EM distributes reads across them and produces a confident-looking split
that is not real. Report group totals, or work at genus level.

### Implausible species usually mean a database gap

An oral epibiont in a potato pit, or a fish pathogen in silage, is the
classifier finding the nearest available reference for something absent from
the database. GTDB resolved both in our data — one became an environmental
*Saccharimonas*, the other an unnamed *Flavobacterium* sp.

---

## Running the pipeline

### The run dies when your SSH session drops

```
[SIGHUP handler] WARN  Killing running tasks (1)
... abortedCount=1
```

Nextflow runs in the foreground on the login node. Closing your laptop kills it.
`-resume` recovers completed tasks, but the in-flight stage is lost — an hour of
basecalling included.

Use `tmux`:

```bash
tmux new -s ont16s
# launch the pipeline, then Ctrl+B then D to detach
tmux attach -t ont16s
```

`nohup ... &` protects against SIGHUP but not against suspension, and a
backgrounded process that tries to read the terminal gets stopped by SIGTTIN.
tmux avoids the whole class of problem.

### Never Ctrl+Z a running pipeline

It suspends the process and leaves the session lock held:

```
ERROR ~ Unable to acquire lock on session with ID ...
```

Recover with:

```bash
jobs
kill -9 %1                              # kill by job number
pkill -9 -u $USER -f nextflow           # or by pattern
lsof <path-to>/db/LOCK                  # no output = stale
rm -f <path-to>/db/LOCK
squeue -M wice -u $USER                 # scancel orphaned SLURM jobs
```

`kill` on a *stopped* process may not land — SIGTERM cannot be delivered until
it resumes — so use `kill -9`. Ctrl+C once is the correct way to stop a run.

### Java disappears between logins

```
/usr/bin/which: no java in (...)
NOTE: Nextflow needs a Java virtual machine to run.
```

Module loads do not persist. `module load Java/...` each session, or install
Nextflow into a conda environment that bundles its own JRE.

---

## Transferring data

### scp cannot resume

A dropped transfer at 80% means restarting from zero. Use rsync:

```bash
caffeinate -i rsync -rvh --progress --partial --append-verify --timeout=300 \
  -e "ssh -o Compression=no -o ServerAliveInterval=20 -o ServerAliveCountMax=30" \
  pod5/ user@cluster:/project/run/pod5/
```

- `--partial --append-verify` resumes a part-transferred file rather than
  restarting it
- `--no-times` if the destination is scratch: preserved mtimes make fresh
  uploads look stale to a 30-day purge policy
- `-o Compression=no` — POD5 is already compressed; SSH compression wastes CPU
- `caffeinate -i` (macOS) stops the machine sleeping mid-transfer

### Diagnose slow transfers before tuning flags

Measure the achievable rate first:

```bash
# download
ssh user@cluster 'dd if=/dev/zero bs=1M count=200 2>/dev/null' | dd of=/dev/null
# upload
dd if=/dev/zero bs=1M count=200 2>/dev/null | ssh user@cluster 'dd of=/dev/null'
```

Links are often asymmetric, so a fine download figure tells you little. If
upload is ~1 MB/s, no rsync flag will help — switch to wired ethernet, drop
VPN, or ask about Globus. At 8 MB/s, ~90 GB takes about 3 hours.

Parallel single-file transfers sometimes beat one stream when per-connection
shaping is the limit:

```bash
ls pod5/*.pod5 | xargs -n 1 -P 4 -I {} \
  rsync -h --partial --append-verify {} user@cluster:/project/run/pod5/
```

### Check for truncated POD5 before spending GPU time

Compare file sizes against timestamps — a file written when the transfer died
will be conspicuously small:

```bash
ls -lh pod5/
for f in pod5/*.pod5; do
  n=$(pod5 view "$f" --include read_id 2>/dev/null | tail -n +2 | wc -l)
  printf "%s\t%s reads\n" "$f" "$n"
done
```

Anything that errors or returns 0 is corrupt. Move it aside; a truncated POD5
fails partway through basecalling, after the GPU allocation is already spent.

### POD5 files are written in time order

Files `_0` through `_5` are all from the first hours of a run. Pore performance
and read-length distribution change over a flow cell's life, so an early-only
subset is not a random sample of the library. Fine for relative abundance,
worth noting as a caveat.
