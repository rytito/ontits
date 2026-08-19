#!/usr/bin/env bash
#
# Construye los entornos conda que el pipeline necesita.
#
#   bash bin/setup_environments.sh /ruta/a/conda_envs
#
# Colócalos en almacenamiento del GRUPO para que todo el equipo comparta una
# copia. Evita el scratch de Lustre: los entornos conda crean ~100k archivos
# pequeños, que Lustre maneja mal, y scratch purga los archivos sin acceso a
# los 30 días.
set -euo pipefail

ENVROOT="${1:?uso: setup_environments.sh <raiz_de_entornos_conda>}"
mkdir -p "$ENVROOT"

MAMBA="${MAMBA:-micromamba}"
command -v "$MAMBA" >/dev/null || {
    echo "ERROR: $MAMBA no se encuentra en el PATH." >&2
    echo "Nota: en muchos clústeres micromamba es una FUNCIÓN del shell, no" >&2
    echo "un ejecutable. Si es el caso, apunta MAMBA= al binario real." >&2
    exit 1
}

echo "==> Entorno ont (manejo de lecturas, demux, QC)"
# cutadapt vive aquí: el demultiplexado dual-index lo necesita, a diferencia
# del pipeline 16S que usaba dorado demux.
"$MAMBA" create -y -p "$ENVROOT/ont" -c conda-forge -c bioconda \
    python=3.11 samtools seqkit minimap2 chopper cutadapt csvtk

# pod5 y NanoPlot van por pip: la receta de bioconda para pod5 fija
# polars >=0.19,<1.dev0, que ya no resuelve contra Python moderno y rompe
# todo el solve.
"$MAMBA" run -p "$ENVROOT/ont" pip install --quiet pod5 NanoPlot || \
    echo "    (pod5/NanoPlot opcionales -- omitidos)"

echo "==> Entorno emu (perfilado taxonómico)"
# El README de Emu apunta a Python 3.7; se instala bajo versiones más
# nuevas, pero fijar 3.10 evita fallos de import ocasionales.
"$MAMBA" create -y -p "$ENVROOT/emu" -c conda-forge -c bioconda \
    python=3.10 emu
"$MAMBA" run -p "$ENVROOT/emu" pip install --quiet osfclient

echo "==> Entorno r (phyloseq)"
"$MAMBA" create -y -p "$ENVROOT/r" -c conda-forge -c bioconda \
    r-base r-tidyverse r-vegan bioconductor-phyloseq

echo
echo "==> Verificación:"
for e in ont emu r; do
    printf "    %-5s " "$e"
    ls "$ENVROOT/$e/bin" 2>/dev/null | head -1 >/dev/null && echo "ok" || echo "FALTA"
done
"$MAMBA" run -p "$ENVROOT/ont" cutadapt --version >/dev/null 2>&1 \
    && echo "    cutadapt ok" || echo "    cutadapt FALTA"
echo
echo "==> Define en conf/site.yml:"
echo "    conda_ont: $ENVROOT/ont"
echo "    conda_emu: $ENVROOT/emu"
echo "    conda_r:   $ENVROOT/r"
