# Uso

> 🌐 **Español** (versión principal) · [English](USAGE.en.md)

La configuración inicial del sitio (Java, Nextflow, Dorado, entorno del
shell, planificación de almacenamiento) es idéntica a la de ont16s — ver
[USAGE de ont16s](../../ont_16S_repo/ont16s/docs/USAGE.md) para el detalle
completo. Aquí, lo específico de ontits.

## Entornos y bases de datos

```bash
bash bin/setup_environments.sh /shared/conda_envs      # añade cutadapt al entorno ont
micromamba activate /shared/conda_envs/emu
bash bin/setup_databases.sh /shared/databases unite-fungi   # o: unite-all | all
```

Las bases preconstruidas de Emu para ITS vienen de UNITE:

| Base de datos | Úsala para | No la uses para |
|---|---|---|
| **unite-fungi** | Taxonomía principal (solo Fungi) | Detectar co-amplificación de plantas |
| **unite-all** | Cuantificar co-amplificación (todos los eucariotas) | — |

Correr ambas es la versión ITS del chequeo GTDB/SILVA de ont16s: KYO2
co-amplifica ITS de plantas y otros eucariotas; contra `unite-fungi` esas
lecturas se fuerzan al hongo más parecido o se pierden. Ambas bases son
pequeñas (≤100 MB): correr las dos cuesta poco.

Las preconstruidas usan el release de UNITE de ~2021. Para la versión
actual, construye con `emu build-database` desde el release general de UNITE.

## Parámetros

### Entrada

| Parámetro | Por defecto | Notas |
|---|---|---|
| `--step` | `all` | `all` \| `from_bam` \| `from_fastq` |
| `--pod5_dir` | — | Con basecalling desactivado, MinKNOW escribe `pod5/`, no `pod5_pass/` |
| `--calls_bam` | — | Para `--step from_bam` |
| `--fastq_dir` | — | Para `--step from_fastq`; archivos `<SampleID>.fastq` **ya sin barcodes/primers** |
| `--metadata` | — | TSV con `SampleID` (`F1R1`…`F4R4`); habilita phyloseq |
| `--outdir` | `results` | |

### Basecalling

Idéntico a ont16s: `--model sup` (→ `sup@v5.2.0` para R10.4.1),
`--max_reads` para pruebas. `--no-trim` está cableado en el pipeline: el
recorte automático de Dorado se comería los barcodes que cutadapt necesita.

### Demultiplexado (cutadapt dual-index)

| Parámetro | Por defecto | Notas |
|---|---|---|
| `--barcodes` | `assets/barcodes.tsv` | Hoja de 16 combinaciones (SampleID, fwd_id, rev_id, fwd_barcode, rev_barcode) |
| `--fwd_primer` | `TAGAGGAAGTAAAAGTCGTAA` | ITS1-F_KYO2 |
| `--rev_primer` | `TTYRCTRCGTTCTTCATC` | ITS2_KYO2 (IUPAC: Y, R) |
| `--demux_error_rate` | `0.15` | ~4 errores tolerados en el adaptador barcode+primer de 30–33 nt |
| `--demux_min_overlap` | `25` | Exige ver barcode **y** primer, no un fragmento |

Dos rondas: (1) `--revcomp` orienta la lectura y asigna el barcode F
(adyacente a su primer, buscados como un único adaptador — con los mismos 4
Golay en ambos extremos, el barcode solo no identifica el extremo); (2)
asigna el barcode R por el 3'. Solo se aceptan lecturas con **ambos**
barcodes; combinaciones fuera de la hoja se descartan como index hopping y
se contabilizan como `invalid_*` en `demux_stats.tsv`.

Si `unassigned_*` supera ~25%: revisa la calidad de la corrida primero;
después baja `--demux_min_overlap` (p. ej. a 20) antes que subir la tasa de
error — los Golay difieren en ≥5 posiciones, pero el margen se come rápido.

### Filtrado y profundidad

| Parámetro | Por defecto | Notas |
|---|---|---|
| `--min_qscore` | `10` | |
| `--min_len` / `--max_len` | `150` / `600` | El amplicón KYO2 mide ~300 pb (inserto ~260 pb sin primers/barcodes); la ventana deja margen para la variación real de ITS1 — no la estreches sin mirar `qc/read_stats.tsv` |
| `--subsample` | `20000` | `0` lo desactiva |
| `--seed` | `42` | Submuestreo reproducible |
| `--min_reads` | `1000` | Más bajo que en 16S; el control negativo (`F4R4`) debe caer por debajo |

**A diferencia del 16S**, la longitud de ITS1 varía biológicamente entre
taxones. Una ventana estrecha sesga la composición contra taxones de ITS
largo o corto — por eso el rango es deliberadamente amplio y solo elimina
fragmentos y concatémeros.

### Taxonomía

| Parámetro | Notas |
|---|---|
| `--emu_dbs` | `nombre:ruta[,...]` — p. ej. `'unite:/db/emu_unite_fungi'` |

## Perfiles

Idénticos a ont16s: `standard`, `slurm`, `vsc_wice`, `vsc_genius`, `conda`,
`apptainer`, `test`, `debug`. Se combinan con comas: `-profile vsc_wice,conda`.

## Control de calidad de la corrida

1. **`results/02_demux/demux_stats.tsv`** — primero. Lecturas por muestra,
   `unassigned_*` (sin uno de los dos barcodes) e `invalid_*` (combinación
   no usada = index hopping medido).
2. **Comunidad mock (`F4R3`)** — composición esperada vs. observada valida
   primers, demux y base de datos de una vez.
3. **Control negativo (`F4R4`)** — debe quedar bajo `--min_reads` y
   excluido del objeto phyloseq. Si no, contaminación.
4. **`07_phyloseq/<db>/uncharacterised_fraction.tsv`** — fracción
   "unidentified"/SH de UNITE por muestra. En sistemas poco estudiados esta
   fracción es un resultado, no un error.

## Ejecución segura y reanudación

Igual que ont16s: corre dentro de `tmux`, nunca Ctrl+Z, reanuda con
`-resume`. Tras un basecalling exitoso y un fallo posterior:

```bash
nextflow run main.nf ... --step from_bam --calls_bam results/01_basecall/calls.bam -resume
```
