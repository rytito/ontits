#!/usr/bin/env bash
#
# Descarga las bases de datos ITS preconstruidas de Emu (UNITE).
# Ejecutar UNA VEZ en un nodo de login -- los nodos de cómputo a menudo no
# tienen internet.
#
#   bash bin/setup_databases.sh /ruta/a/databases [unite-fungi|unite-all|all]
#
set -euo pipefail

DBROOT="${1:?uso: setup_databases.sh <raiz_de_bases> [cual]}"
WHICH="${2:-unite-fungi}"

mkdir -p "$DBROOT"
echo "==> Raíz de bases de datos: $DBROOT"
echo "    Usa almacenamiento de PROYECTO, no scratch: scratch purga archivos"
echo "    sin acceso por 30 días, y estas son costosas de re-descargar."
echo

command -v osf >/dev/null || pip install --quiet osfclient
echo "    Bases de datos preconstruidas disponibles en el OSF de Emu:"
osf -p 56uf7 ls | grep -i 'prebuilt.*unite' || true
echo

# ------------------------------------------------------------ UNITE fungi
# UNITE restringida a Fungi: la referencia estándar para metabarcoding ITS.
# Es la base principal para este pipeline. Pequeña (~50 MB), Emu corre
# holgado en 16-32 GB de RAM.
if [[ "$WHICH" == "unite-fungi" || "$WHICH" == "all" ]] && [ ! -d "$DBROOT/emu_unite_fungi" ]; then
    echo "==> UNITE fungi"
    mkdir -p "$DBROOT/emu_unite_fungi" && cd "$DBROOT/emu_unite_fungi"
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/unite-fungi.tar
    tar -xf unite-fungi.tar && rm unite-fungi.tar
fi

# ------------------------------------------------------------ UNITE all
# UNITE con todos los eucariotas. Útil como control de contaminación: el
# par KYO2 co-amplifica ITS de plantas y otros eucariotas, y contra
# unite-fungi esas lecturas se forzarían hacia el hongo más parecido o se
# perderían. Correr ambas y comparar es la versión ITS del chequeo
# GTDB/SILVA del pipeline 16S.
if [[ "$WHICH" == "unite-all" || "$WHICH" == "all" ]] && [ ! -d "$DBROOT/emu_unite_all" ]; then
    echo "==> UNITE all (todos los eucariotas)"
    mkdir -p "$DBROOT/emu_unite_all" && cd "$DBROOT/emu_unite_all"
    osf -p 56uf7 fetch osfstorage/emu-prebuilt/unite-all.tar
    tar -xf unite-all.tar && rm unite-all.tar
fi

echo
echo "==> Listo. Verifica que cada base tenga species_taxid.fasta y taxonomy.tsv:"
find "$DBROOT" -maxdepth 2 -name "taxonomy.tsv" | sed 's/^/    /'
echo
echo "==> Luego define en conf/site.yml:"
echo "    emu_dbs: 'unite:$DBROOT/emu_unite_fungi'"
echo "    # o ambas:"
echo "    # emu_dbs: 'unite:$DBROOT/emu_unite_fungi,unite_all:$DBROOT/emu_unite_all'"
echo
echo "NOTA: las preconstruidas de Emu usan UNITE 2021 aprox. Si necesitas la"
echo "versión actual de UNITE, construye la base con 'emu build-database' a"
echo "partir del release general FASTA de UNITE y su taxonomía."
