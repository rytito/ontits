# Publishing to GitHub

> 🌐 [Español](SETUP_GITHUB.md) (primary version) · **English**

## 1. Create the repository locally

```bash
cd ontits
git init
git add .
git commit -m "ontits v1.0.0: ONT ITS dual-index pipeline, POD5 to phyloseq"
git branch -M main
```

`conf/site.yml` and `metadata.tsv` are gitignored — only `site.yml.example`
and `assets/metadata_template.tsv` are versioned, since the real files carry
absolute paths, an account name and run-specific sample metadata. The
`assets/barcodes.tsv` sheet IS versioned: it is part of the library design,
not a secret.

**Do not push the parent folder** (`ont_ITS/`): it holds lab worksheets
(`primers.xlsx`) that do not belong in a public repository.

## 2. Push

```bash
gh repo create ontits --public --source=. --push
```

Or after creating it in the web UI:

```bash
git remote add origin git@github.com:<user>/ontits.git
git push -u origin main
git tag -a v1.0.0 -m "First release, demux validated on synthetic reads"
git push --tags
```

## 3. Deploy to a cluster as a clone, not a copy

If the cluster copy and the repository are independent they diverge, and you
end up not knowing which one is authoritative.

```bash
cd /shared/pipelines
git clone https://github.com/<user>/ontits.git
cd ontits
cp conf/site.yml.example conf/site.yml
$EDITOR conf/site.yml
```

Updates then arrive with `git pull`, and nobody edits files in place.

Or let Nextflow handle it — it clones and caches automatically:

```bash
nextflow run <user>/ontits -r v1.0.0 \
    -profile vsc_wice,conda -params-file /shared/pipelines/site.yml \
    --step all --pod5_dir /project/run/pod5 --metadata metadata.tsv
```

Pinning `-r v1.0.0` keeps a colleague's run reproducible even after you push
changes. Note: with `nextflow run <user>/ontits`, the `site.yml` and
`metadata.tsv` live OUTSIDE the cached clone — keep them at a stable path and
pass them as absolute paths.

## 4. Smoke test before telling anyone

Cheapest first — compiles the script, submits nothing:

```bash
nextflow run main.nf --help
```

Then the full path with `-profile test`, which caps basecalling at 20k reads:

```bash
nextflow run main.nf \
    -profile vsc_wice,conda,test -params-file conf/site.yml \
    --step all --pod5_dir /project/run/pod5 --metadata metadata.tsv \
    --outdir results_test
```

When it finishes, check `results_test/02_demux/demux_stats.tsv`: most reads
should land on the 16 samples, not in `unassigned_*`.

## 5. What to tell your team

> Pipeline: `github.com/<user>/ontits`, deployed at `/shared/pipelines/ontits`.
>
> `module load Java/...`, then run inside `tmux` — otherwise Nextflow dies
> with your SSH session.
>
> **Read `docs/GOTCHAS.md` first.** Several documented failures report
> success while silently losing data.
>
> Keep POD5 on project storage, and the UNITE databases on data/project —
> scratch purges after 30 days.
