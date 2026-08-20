# ontits

> 🌐 **Español** (versión principal) · [English](README.en.md)

Pipeline de Nextflow para secuenciación Oxford Nanopore de amplicones **ITS
fúngicos con dual-index estilo Illumina**: desde la señal cruda en POD5 hasta
un objeto `phyloseq`, en un solo comando.

Derivado de [ont16s](../ont_16S_repo/) y construido para la **química Kit 14**
(R10.4.1 / FLO-MIN114) con el basecalling **desactivado** en el secuenciador,
ejecutado en un clúster HPC.

```
POD5  ──▶  Dorado basecall      (GPU, sup@v5.2.0, --no-trim)
      ──▶  samtools fastq       (BAM completo -> un FASTQ)
      ──▶  cutadapt demux       (dual-index Golay 4F x 4R = 16 muestras,
                                 ambos barcodes obligatorios, recorta
                                 adaptadores + barcodes + primers)
      ──▶  chopper              (q≥10, 150–600 pb; amplicón ~300 pb)
      ──▶  seqkit sample        (normalización de profundidad)
      ──▶  Emu abundance        (UNITE fungi; opcional UNITE all)
      ──▶  Emu combine          (tablas taxones × muestras)
      ──▶  phyloseq             (.rds + tablas de QC + gráficos)
```

## El esquema de amplicón

La librería **no** usa un kit de barcoding de ONT. Usa los primers de Toju
*et al.* 2012 con construcciones dual-index tipo Earth Microbiome, secuenciadas
en Nanopore:

```
5'─[P5 + sitio Read1 (58 nt)]─[Golay F (12 nt)]─[ITS1-F_KYO2 (21 nt)]─
      ...inserto ITS...
   ─[rc ITS2_KYO2 (18 nt)]─[rc Golay R (12 nt)]─[rc sitio Read2 + P7 (58 nt)]─3'
```

- **Primer F:** ITS1-F_KYO2 `TAGAGGAAGTAAAAGTCGTAA`
- **Primer R:** ITS2_KYO2 `TTYRCTRCGTTCTTCATC`
- **Barcodes Golay (12 nt), los mismos 4 en ambos extremos:**
  `AGCCTTCGTCGC` (1), `TCCATACCGGAA` (2), `AGCCCTGCTACA` (3), `CCTAACGGTCCA` (4)
- **16 combinaciones** = 16 muestras, definidas en
  [`assets/barcodes.tsv`](assets/barcodes.tsv). En la corrida original,
  `F4R3` es la comunidad mock y `F4R4` el control negativo.

Por eso `dorado demux` no aplica: el demultiplexado corre con **cutadapt** en
dos rondas. La ronda 1 orienta cada lectura (`--revcomp`) y la asigna por
`barcodeF+primerF` recortando el 5'; la ronda 2 la asigna por
`rc(barcodeR+primerR)` recortando el 3'. Como los mismos 4 Golay van en ambos
extremos, un barcode aislado no identifica el extremo — el primer adyacente sí,
así que barcode y primer se buscan como un único adaptador. Solo se aceptan
lecturas con **ambos** barcodes (el equivalente de `--barcode-both-ends`);
las combinaciones fuera de la hoja se descartan como index hopping.

## Por qué Emu + UNITE

- **Emu** resuelve lecturas ambiguas mediante expectation-maximisation sobre
  alineamientos de longitud completa — el mismo razonamiento que en 16S aplica
  a ITS, donde las bases de referencia están llenas de secuencias casi
  idénticas entre especies hermanas.
- **UNITE** es la referencia estándar para metabarcoding ITS. Emu distribuye
  dos preconstruidas: `unite-fungi` (solo hongos — la base principal) y
  `unite-all` (todos los eucariotas — control de co-amplificación: KYO2
  también amplifica ITS de plantas, y contra `unite-fungi` esas lecturas se
  forzarían al hongo más parecido).
- **TSV planos y un objeto phyloseq**, no un reporte del que haya que extraer
  los datos a mano.

## Inicio rápido

```bash
# 1. Entornos y bases de datos (una vez por sitio)
bash bin/setup_environments.sh /shared/conda_envs
micromamba activate /shared/conda_envs/emu
bash bin/setup_databases.sh /shared/databases unite-fungi

# 2. Configuración del sitio
cp conf/site.yml.example conf/site.yml
$EDITOR conf/site.yml          # rutas, cuenta, bases de datos

# 3. Ejecutar
nextflow run main.nf \
    -profile slurm,conda \
    -params-file conf/site.yml \
    --step all \
    --pod5_dir /project/myrun/pod5 \
    --metadata metadata.tsv \
    --outdir results
```

**Ejecútalo dentro de `tmux`.** Nextflow corre en primer plano y muere junto
con tu sesión SSH — `-resume` recupera lo completado, pero se pierde la etapa
en curso.

## Etapas

Se seleccionan con `--step`. La sintaxis estricta de Nextflow 26.04 eliminó
`-entry`.

| Punto de partida | Comando |
|---|---|
| POD5 (corrida completa) | `--step all --pod5_dir <dir>` |
| Un `calls.bam` existente | `--step from_bam --calls_bam calls.bam` |
| FASTQ por muestra, ya recortados | `--step from_fastq --fastq_dir <dir>` |

`from_bam` es lo que necesitas tras un fallo de demultiplexado o de perfilado —
omite por completo la etapa costosa de GPU.

## Salidas

```
results/
├── 01_basecall/        versiones de las herramientas
├── 02_demux/           demux_stats.tsv (lecturas por muestra, no asignadas,
│                       combinaciones inválidas) + logs de cutadapt
├── 03_filtered/        FASTQ por muestra, sin barcodes/primers, filtrado
├── 04_subsampled/      FASTQ con profundidad normalizada
├── 05_emu/<db>/        tablas de abundancia por muestra
├── 06_tables/          <db>_emu-combined-{species,genus,family}[-counts].tsv
├── 07_phyloseq/<db>/   phyloseq.rds, diversidad alfa, ordenación, gráficos
├── qc/                 estadísticas de lecturas, lecturas por muestra
└── pipeline_info/      timeline, reporte, trace, DAG
```

Revisa `02_demux/demux_stats.tsv` antes de nada: si una fracción grande cae en
`unassigned_*`, baja `--demux_min_overlap` o sube `--demux_error_rate` con
cautela — y si aparecen lecturas en `invalid_*`, eso es index hopping real.

El objeto se carga directamente:

```r
ps <- readRDS("results/07_phyloseq/unite/phyloseq.rds")
```

## Metadatos

Un TSV con una columna `SampleID` que coincida con los nombres de la hoja de
barcodes (`F1R1` … `F4R4`). Ver `assets/metadata_template.tsv`. Cualquier otra
columna se convierte en variable de muestra, y las categóricas se prueban
automáticamente con PERMANOVA.

```tsv
SampleID	sample_name	group	replicate
F1R1	S01	group_a	1
F4R3	MockCom	control	1
```

Sin `--metadata` el pipeline se detiene después de `06_tables`.

**Usa la comunidad mock (`F4R3`)** para validar la corrida: es el
ZymoBIOMICS Microbial Community Standard (Zymo D6300), y con primers ITS
solo sus dos levaduras amplifican — el resultado esperado es *Saccharomyces
cerevisiae* + *Cryptococcus neoformans* y nada más (ver docs/USAGE.md).
**El blanco (`F4R4`)** debe caer por debajo de `--min_reads` y quedar
excluido; si no, investiga contaminación antes de interpretar nada.

## Requisitos

- Nextflow **≥ 26.04** y Java 11+ (`module load Java/...` en la mayoría de
  clústeres)
- Una GPU NVIDIA con tensor cores. **Dorado no soporta P100/GP100.**
- Dorado ≥ 2.0 (binario estático; deliberadamente sin contenedor)
- conda/micromamba, o Apptainer
- cutadapt ≥ 4 (lo instala `setup_environments.sh` en el entorno `ont`)

## Documentación

- [`docs/USAGE.md`](docs/USAGE.md) — parámetros, perfiles, configuración del
  clúster
- [`docs/GOTCHAS.md`](docs/GOTCHAS.md) — modos de fallo que encontramos, y cómo
  detectarlos

**Lee `GOTCHAS.md` antes de tu primera corrida real.** Varios de los fallos
documentados reportan éxito mientras pierden datos silenciosamente.

## Cómo citar

Cita las herramientas subyacentes: Dorado (Oxford Nanopore), cutadapt (Martin
2011), Emu (Curry *et al.* 2022, *Nature Methods*), minimap2 (Li 2018),
chopper, seqkit, phyloseq (McMurdie & Holmes 2013), los primers ITS1-F_KYO2 /
ITS2_KYO2 (Toju *et al.* 2012, *PLoS ONE*), más **UNITE con su número de
versión** (Nilsson *et al.* 2019; las preconstruidas de Emu usan el release
2021 aprox. — verifica y repórtalo en tus métodos).

Si reportas resultados de UNITE, indica que los placeholders tipo
"unidentified" y las hipótesis de especie (SH) son parte del diseño de la base
de datos: en un sistema poco estudiado esa fracción es un hallazgo, no un
error.

## Financiamiento y créditos

Este pipeline fue desarrollado dentro del proyecto **MicroAndes** —
*Aprovechar los microbios para fortalecer a las comunidades andinas de
microproductores de alimentos fermentados* — una **Short Initiative de
VLIR-UOS** (proyecto KU Leuven
[3M250529](https://research.kuleuven.be/portal/nl/project/3M250529),
referencia de financiamiento PE2025SIN468A101, 2025–2027).

- **Institución coordinadora:** KU Leuven — Laboratorio de Bacteriología
  Molecular, Instituto Rega
- **Instituciones socias:** Universidad Nacional del Centro del Perú (UNCP) y
  comunidades andinas de microproductores
- **Financiamiento:** VLIR-UOS (Vlaamse Interuniversitaire Raad), programa
  Short Initiatives

## Licencia

MIT — ver [`LICENSE`](LICENSE).
