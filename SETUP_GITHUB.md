# Publicación en GitHub

> 🌐 **Español** (versión principal) · [English](SETUP_GITHUB.en.md)

## 1. Crear el repositorio localmente

```bash
cd ontits
git init
git add .
git commit -m "ontits v1.0.0: ONT ITS dual-index pipeline, POD5 to phyloseq"
git branch -M main
```

`conf/site.yml` y `metadata.tsv` están en el gitignore — solo se versionan
`site.yml.example` y `assets/metadata_template.tsv`, ya que los archivos
reales contienen rutas absolutas, un nombre de cuenta y metadatos de muestras
específicos de una corrida. La hoja `assets/barcodes.tsv` SÍ se versiona: es
parte del diseño de la librería, no un secreto.

**No subas la carpeta madre** (`ont_ITS/`): contiene hojas de laboratorio
(`primers.xlsx`) que no pertenecen a un repositorio público.

## 2. Subir

```bash
gh repo create ontits --public --source=. --push
```

O después de crearlo en la interfaz web:

```bash
git remote add origin git@github.com:<user>/ontits.git
git push -u origin main
git tag -a v1.0.0 -m "First release, demux validated on synthetic reads"
git push --tags
```

## 3. Desplegar en un clúster como clon, no como copia

Si la copia del clúster y el repositorio son independientes, divergen y
terminas sin saber cuál es la autoritativa.

```bash
cd /shared/pipelines
git clone https://github.com/<user>/ontits.git
cd ontits
cp conf/site.yml.example conf/site.yml
$EDITOR conf/site.yml
```

Las actualizaciones llegan entonces con `git pull`, y nadie edita archivos en
el lugar.

O deja que Nextflow lo maneje — clona y cachea automáticamente:

```bash
nextflow run <user>/ontits -r v1.0.0 \
    -profile vsc_wice,conda -params-file /shared/pipelines/site.yml \
    --step all --pod5_dir /project/run/pod5 --metadata metadata.tsv
```

Fijar `-r v1.0.0` hace que la corrida de un colega sea reproducible incluso
después de que subas cambios. Nota: con `nextflow run <user>/ontits`, el
`site.yml` y el `metadata.tsv` viven FUERA del clon cacheado — guárdalos en
una ruta estable y pásalos con rutas absolutas.

## 4. Prueba de humo antes de avisar a nadie

Lo más barato primero — compila el script, no envía nada:

```bash
nextflow run main.nf --help
```

Luego la ruta completa con `-profile test`, que limita el basecalling a 20k
lecturas:

```bash
nextflow run main.nf \
    -profile vsc_wice,conda,test -params-file conf/site.yml \
    --step all --pod5_dir /project/run/pod5 --metadata metadata.tsv \
    --outdir results_test
```

Al terminar, revisa `results_test/02_demux/demux_stats.tsv`: la mayoría de
las lecturas debe caer en las 16 muestras, no en `unassigned_*`.

## 5. Qué decirle a tu equipo

> Pipeline: `github.com/<user>/ontits`, desplegado en `/shared/pipelines/ontits`.
>
> `module load Java/...`, y luego correr dentro de `tmux` — de lo contrario
> Nextflow muere con tu sesión SSH.
>
> **Lean `docs/GOTCHAS.md` primero.** Varios de los fallos documentados
> reportan éxito mientras pierden datos silenciosamente.
>
> Mantengan los POD5 en el almacenamiento del proyecto, y las bases UNITE en
> data/proyecto — scratch purga a los 30 días.
