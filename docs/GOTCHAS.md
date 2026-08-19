# Gotchas (trampas conocidas)

> **Nota (ontits):** heredado de ont16s. Todo lo relativo a cuotas de disco, Dorado, SLURM/VSC, Nextflow, Emu, tmux y transferencia de datos aplica tal cual. Las secciones sobre `dorado demux` / `--kit-name`, GTDB/SILVA y resolución de especies por 16S NO aplican aquí: ontits demultiplexa con cutadapt (ver README) y perfila contra UNITE. Las lecciones sobre combinaciones de barcodes sin usar como medición de cross-talk y sobre muestras de baja profundidad aplican doblemente a ITS.

> 🌐 **Español** (versión principal) · [English](GOTCHAS.en.md)

Modos de fallo que realmente encontramos, con la forma de detectar cada uno.
Varios reportan éxito mientras pierden datos silenciosamente.

---

## Pérdida silenciosa de datos

### La cuota de disco trunca los BAM, y Dorado igual sale con código 0

**Síntoma:** `EOF marker is absent. The input is probably truncated`, mientras
`sacct` muestra `COMPLETED 0:0`.

**Causa:** el scratch de VSC es de **100 GB**, no 500. Con ~86 GB de POD5 más
`calls.bam` más la salida del demux, las escrituras chocan contra el límite.
Dorado no propaga el fallo de escritura como código de salida distinto de cero.

**Detección:**
```bash
myquota
samtools quickcheck -u -v calls.bam
find demux -name '*.bam' | xargs samtools quickcheck -u -v && echo "ALL OK"
```

El `-u` importa: sin él, `quickcheck` reporta `had no targets in header`, lo
cual es **normal** en un BAM no alineado. Solo "missing EOF block" significa
truncamiento.

**Solución:** mantén los POD5 en el almacenamiento del proyecto y apunta
`--pod5_dir` allí. El pipeline establece `stageInMode = 'symlink'` para que
Nextflow no los duplique.

### Leer salidas mientras un trabajo todavía escribe

Los mismos errores de truncamiento aparecen si cuentas lecturas en `demux/`
mientras `dorado demux` sigue corriendo. Revisa `squeue -M wice -u $USER` antes
de concluir que algo está roto.

---

## Dorado

### `--kit-name` en el basecaller no separa las lecturas

Solo las *etiqueta*. `dorado demux` sigue siendo necesario para producir
archivos por muestra. El malentendido más común sobre la herramienta.

### `--no-trim` es obligatorio si se demultiplexa después

Sin él, el basecaller recorta adaptadores y barcodes, y el demux no encuentra
nada — todo termina en `unclassified`.

### La salida es un árbol anidado estilo MinKNOW

```
demux/<experiment>/<sample_id>/<run_id>/bam_pass/barcodeNN/*.bam
```

Reconstruido desde los metadatos del POD5. Cuando MinKNOW no tenía sample ID
configurado, el nivel intermedio es literalmente `unknown`. **Nunca uses glob
con una ruta fija** — recorre con `find`. Cualquier script que use
`demux/*.fastq` no coincide con nada, silenciosamente.

### Las GPU P100 no están soportadas

Los requisitos de Dorado excluyen P100/GP100. En Genius esta suele ser la
partición más vacía y por lo tanto la más tentadora. Usa `gpu_v100`, o wICE.

### `ReqNodeNotAvail, Reserved for maintenance` sin reserva alguna

`scontrol show reservation` puede no reportar nada mientras el trabajo sigue
sin programarse — los nodos están en drenado.

```bash
sinfo -M wice -p gpu_h100 -o "%P %a %D %t %N"
sinfo -M wice -R
```

`maint` y un `$` al final significan fuera de servicio; `mix`/`alloc` solo
significan ocupado. `SCHEDNODES (null)` con `START_TIME N/A` significa que
ningún nodo podría satisfacer jamás la solicitud.

---

## SLURM en VSC

### `--nodes` debe ser explícito en los envíos con GPU

```
sbatch: error: [ERROR] Please specify the number of nodes (using --nodes)
```

Nextflow deriva `--nodes` de `cpus` para los procesos ordinarios, pero un
override de `clusterOptions` reemplaza la cadena completa. La etiqueta de GPU
debe llevar `--nodes=1 --ntasks=1 --gpus-per-node=1` explícitamente.

### Techos de cores por GPU

Excederlos hace que el trabajo sea rechazado de plano.

| Clúster | Partición | Máx. cores/GPU | Máx. mem/GPU |
|---|---|---|---|
| wICE | `gpu_h100` | 16 | 187200 MiB |
| wICE | `gpu_a100` | 18 | 126000 MiB |
| Genius | `gpu_v100` | 4 | 84000 MiB |

### `-M` / `--clusters` es obligatorio

Un `squeue` a secas apunta por defecto al clúster del nodo de login y se ve
vacío incluso con trabajos corriendo. Siempre `squeue -M wice -u $USER`.

### `--parsable` devuelve `jobid;cluster`

Por eso `--dependency=afterok:$JID` falla con "Job dependency problem":

```bash
RAW=$(sbatch --parsable job.slurm)
JID=${RAW%%;*}
```

### `#!/bin/bash -l` es obligatorio

Sin `-l` el shell de login nunca se inicializa, el módulo del clúster no se
carga, y las adiciones al `PATH` de `~/.bashrc` están ausentes.

### Un `srun` en segundo plano no sobrevive al logout

Una asignación suspendida con Ctrl+Z y reanudada con `bg` muere cuando la
sesión SSH se corta. Usa `sbatch`.

### Pegar un bloque mientras `srun` está en cola lo ejecuta en el nodo de login

El `srun` espera; cada línea siguiente queda en el buffer del shell de login y
se ejecuta allí. Así es como Dorado termina reportando "no CUDA devices
available" en un nodo de login. Envía una línea, espera a que cambie el prompt,
confirma con `hostname` y `nvidia-smi`, y luego continúa.

---

## Nextflow

### `-entry` eliminado en 26.04

```
ERROR ~ The `-entry` option is not supported with the strict parser
```

La sintaxis estricta (parser v2) pasó a ser la predeterminada. Este pipeline
usa `--step` en su lugar — un parámetro del pipeline, con doble guion:

```bash
nextflow run main.nf --step from_fastq --fastq_dir in
```

### Otros rechazos de la sintaxis estricta

Los cuatro aparecieron durante la migración:

| Construcción | Reemplazo |
|---|---|
| `exit 1, "msg"` | `error "msg"` |
| `switch`/`case` | cadena de `if` / `else if` |
| `workflow.onComplete` | elimínalo; Nextflow imprime su propio resumen |
| sentencias de nivel superior | muévelas dentro de un workflow o función |

El `exit 1` de Bash **dentro** de un bloque `script:` no se ve afectado — solo
el código Groovy lo está.

### `publishDir` con variables de entrada necesita un closure

```groovy
publishDir "${params.outdir}/${db_name}"      // No such variable: db_name
publishDir { "${params.outdir}/${db_name}" }  // correcto
```

Las directivas se evalúan antes de que las variables de entrada se enlacen. Lo
mismo aplica a `saveAs`, donde el parámetro del closure debe nombrarse
explícitamente (`{ fn -> ... }`) en lugar de confiar en el `it` implícito.

### Los parámetros de línea de comandos llegan como Strings

```
Cannot compare java.lang.String with value '5000' and java.lang.Integer
```

`--subsample 5000` es un String; `params.subsample > 0` entonces falla. Usa
`params.subsample.toString().toInteger()`. Añadir `nf-schema` con declaraciones
tipadas lo manejaría correctamente.

### Nunca definas `PATH` en un bloque `env {}`

Nextflow emite los exports de `env` **después** de `beforeScript`, así que una
asignación de `PATH` allí sobrescribe silenciosamente lo que `beforeScript`
haya definido. Síntoma: `.command.run` contiene dos líneas `export PATH=`, y la
segunda no incluye tus herramientas. Pon las rutas de herramientas solo en
`process.beforeScript`.

### conda y micromamba son *funciones* del shell en VSC

```bash
which micromamba   # imprime el cuerpo de una función, no una ruta
```

Los procesos de Nextflow nunca cargan `~/.bashrc`, así que la directiva `conda`
no puede encontrarlos. Pon el `bin/` de cada entorno en el `PATH` vía
`beforeScript` en su lugar.

### `withName: '!EMU_.*'` no es válido

Los selectores de nombre negados no coinciden con nada, silenciosamente. Define
un valor por defecto a nivel de proceso y sobreescribe para el caso específico
— los selectores posteriores y más específicos ganan.

### No hagas Ctrl+Z a un pipeline en ejecución

Deja retenido el lock de sesión:

```
ERROR ~ Unable to acquire lock on session with ID ...
```

Se arregla con `jobs` y luego `kill %1`. Usa Ctrl+C una vez en su lugar —
Nextflow cancela limpiamente y libera el lock.

---

## Emu

### SILVA necesita ~48 GB de RAM

SILVA 138.2 tiene **2.2 millones de secuencias de referencia (3.2 Gpb)**. Solo
el índice de minimap2 supera los 16 GB, y una asignación de 16 GB muere por OOM:

```
Command 'minimap2 -ax map-ont ...' died with <Signals.SIGKILL: 9>
```

La etiqueta `cpu_med` solicita 48 GB y los duplica al reintentar. GTDB (77k
secuencias) corre cómodamente en 16 GB. Incluir SILVA aproximadamente triplica
el costo de una corrida — vale la pena para el chequeo de cloroplastos,
innecesario para el perfilado rutinario.

### Conteos fraccionarios

EM produce estimaciones como `39879.125`. Redondea antes de DESeq2 o de
`rarefy_even_depth()`, que necesitan enteros.

### `--keep-counts` es obligatorio para `combine-outputs --counts`

Sin él solo obtienes abundancias relativas.

### El tiempo de ejecución escala fuertemente con la profundidad

90k lecturas contra SILVA tomaron >20 min/muestra y 284 iteraciones de EM,
frente a 32 contra la base de datos por defecto. Submuestrea a 20k — las
abundancias relativas son estables muy por debajo de eso.

### Versión de Python

El README de Emu apunta a Python 3.7; bioconda lo instala bajo 3.13. Funciona,
pero exporta `PYTHONNOUSERSITE=1` para que los paquetes de `~/.local` no
eclipsen el entorno. Si fallan los imports, reconstruye fijando 3.10.

---

## Bases de datos

### Las especies de GTDB son binomios — no dividas los encabezados por espacios

```
>RS_GCF_000566285.1 d__Bacteria;...;s__Flavobacterium fluviatile [ssu_len=1510]
```

Un `header.split()` ingenuo trunca la especie a `Flavobacterium` y colapsa
todas las especies de un género en un solo taxón. **Síntoma:** muchos menos
taxones únicos que secuencias — nosotros vimos 20,910 de 75,602. Extrae desde
`d__` hasta el primer ` [`. `bin/gtdb_to_emu.py` lo hace correctamente.

**Verifica después de construir:**
```bash
head -2 emu_input/taxonomy.tsv | awk -F'\t' '{print $2}'
# se espera "Flavobacterium fluviatile", no "Flavobacterium"
```

### `column -t` desfigura la tabla de taxonomía

Divide por cualquier espacio en blanco, así que el espacio dentro de un binomio
parece un salto de columna y todos los rangos aparecen desplazados. El archivo
está bien:

```bash
head -3 taxonomy.tsv | awk -F'\t' '{for(i=1;i<=NF;i++) printf "[%s] ", $i; print ""}'
```

### Usa `ssu_reps`, no `ssu_all`

`ssu_all` contiene todas las copias de SSU de 732k genomas. `ssu_reps` es una
por clúster de especie — ~77k secuencias, comparable a la base por defecto de
Emu.

### Los sufijos de género de GTDB son reales

`Clostridium_S`, `Clostridium_B`, `Clostridium_AM` son **géneros distintos**,
añadidos donde el género de NCBI es polifilético. Nunca los elimines con una
regex — fusiona silenciosamente linajes distintos. Menciónalos en tus métodos
para que los revisores no los lean como erratas.

### GTDB es solo de procariotas

Sin referencias de cloroplastos ni mitocondrias. Las lecturas de plantas mapean
a la cianobacteria más cercana y reciben un nombre de apariencia plausible. Usa
SILVA para cuantificar esa fracción.

---

## Interpretación

### Las muestras de baja profundidad muestran la *máxima* diversidad aparente

Una muestra con 669 lecturas dio 27 taxones por encima del 1% — el máximo de la
corrida — mientras que una muestra sana dio 7. Máxima diversidad a mínima
profundidad es la firma del ruido. Descarta muestras por debajo de ~5000
lecturas.

### Los barcodes sin usar son tu única medición de cross-talk

Llenar los 24 espacios te deja sin control negativo y sin piso de detección.
Deja uno vacío, o incluye un blanco de extracción.

### Las especies estrechamente emparentadas no se separan por 16S

El grupo *L. casei* y los clostridios solventogénicos tienen genes 16S casi
idénticos. EM distribuye las lecturas entre ellos y produce una separación de
apariencia confiable que no es real. Reporta totales por grupo, o trabaja a
nivel de género.

### Las especies inverosímiles suelen indicar un vacío en la base de datos

Un epibionte oral en un pozo de papa, o un patógeno de peces en ensilaje, es el
clasificador encontrando la referencia disponible más cercana para algo ausente
de la base de datos. GTDB resolvió ambos en nuestros datos — uno resultó ser un
*Saccharimonas* ambiental, el otro un *Flavobacterium* sp. sin nombre.

---

## Ejecución del pipeline

### La corrida muere cuando se corta tu sesión SSH

```
[SIGHUP handler] WARN  Killing running tasks (1)
... abortedCount=1
```

Nextflow corre en primer plano en el nodo de login. Cerrar tu laptop lo mata.
`-resume` recupera las tareas completadas, pero la etapa en curso se pierde —
incluida una hora de basecalling.

Usa `tmux`:

```bash
tmux new -s ont16s
# lanza el pipeline, luego Ctrl+B y D para desacoplar
tmux attach -t ont16s
```

`nohup ... &` protege contra SIGHUP pero no contra la suspensión, y un proceso
en segundo plano que intenta leer la terminal es detenido por SIGTTIN. tmux
evita toda esta clase de problemas.

### Nunca hagas Ctrl+Z a un pipeline en ejecución

Suspende el proceso y deja retenido el lock de sesión:

```
ERROR ~ Unable to acquire lock on session with ID ...
```

Recupera con:

```bash
jobs
kill -9 %1                              # matar por número de trabajo
pkill -9 -u $USER -f nextflow           # o por patrón
lsof <ruta-a>/db/LOCK                   # sin salida = obsoleto
rm -f <ruta-a>/db/LOCK
squeue -M wice -u $USER                 # scancel a los trabajos SLURM huérfanos
```

`kill` sobre un proceso *detenido* puede no surtir efecto — SIGTERM no puede
entregarse hasta que se reanude — así que usa `kill -9`. Ctrl+C una vez es la
forma correcta de detener una corrida.

### Java desaparece entre sesiones

```
/usr/bin/which: no java in (...)
NOTE: Nextflow needs a Java virtual machine to run.
```

Los módulos cargados no persisten. `module load Java/...` en cada sesión, o
instala Nextflow en un entorno conda que incluya su propio JRE.

---

## Transferencia de datos

### scp no puede reanudar

Una transferencia cortada al 80% significa reiniciar desde cero. Usa rsync:

```bash
caffeinate -i rsync -rvh --progress --partial --append-verify --timeout=300 \
  -e "ssh -o Compression=no -o ServerAliveInterval=20 -o ServerAliveCountMax=30" \
  pod5/ user@cluster:/project/run/pod5/
```

- `--partial --append-verify` reanuda un archivo parcialmente transferido en
  lugar de reiniciarlo
- `--no-times` si el destino es scratch: los mtimes preservados hacen que las
  subidas recientes parezcan antiguas ante una política de purga de 30 días
- `-o Compression=no` — POD5 ya está comprimido; la compresión SSH desperdicia
  CPU
- `caffeinate -i` (macOS) evita que la máquina se duerma a mitad de la
  transferencia

### Diagnostica las transferencias lentas antes de ajustar flags

Mide primero la tasa alcanzable:

```bash
# descarga
ssh user@cluster 'dd if=/dev/zero bs=1M count=200 2>/dev/null' | dd of=/dev/null
# subida
dd if=/dev/zero bs=1M count=200 2>/dev/null | ssh user@cluster 'dd of=/dev/null'
```

Los enlaces suelen ser asimétricos, así que una buena cifra de descarga dice
poco. Si la subida es ~1 MB/s, ningún flag de rsync ayudará — cambia a ethernet
por cable, quita la VPN, o pregunta por Globus. A 8 MB/s, ~90 GB toman unas 3
horas.

Las transferencias paralelas de archivos individuales a veces superan a un solo
flujo cuando el límite es el shaping por conexión:

```bash
ls pod5/*.pod5 | xargs -n 1 -P 4 -I {} \
  rsync -h --partial --append-verify {} user@cluster:/project/run/pod5/
```

### Revisa si hay POD5 truncados antes de gastar tiempo de GPU

Compara los tamaños de archivo contra las marcas de tiempo — un archivo escrito
cuando la transferencia murió será notoriamente pequeño:

```bash
ls -lh pod5/
for f in pod5/*.pod5; do
  n=$(pod5 view "$f" --include read_id 2>/dev/null | tail -n +2 | wc -l)
  printf "%s\t%s reads\n" "$f" "$n"
done
```

Cualquiera que dé error o devuelva 0 está corrupto. Apártalo; un POD5 truncado
falla a mitad del basecalling, cuando la asignación de GPU ya se gastó.

### Los archivos POD5 se escriben en orden temporal

Los archivos `_0` a `_5` son todos de las primeras horas de una corrida. El
rendimiento de los poros y la distribución de longitudes de lectura cambian a
lo largo de la vida de una celda de flujo, así que un subconjunto solo inicial
no es una muestra aleatoria de la librería. Aceptable para abundancia relativa,
vale la pena anotarlo como salvedad.
