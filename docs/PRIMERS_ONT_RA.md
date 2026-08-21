# Rediseño de primers ITS para Rapid Adapter (RA) — SQK-RPB114.24

> 🌐 Español (versión principal) · English summary at the end

## Por qué el rediseño

El diseño original usaba construcciones dual-index estilo Illumina (P5/P7 +
barcodes Golay). Dos problemas para Nanopore con Rapid Adapter:

1. Los handles de Illumina son inertes en ONT: 116 nt de oligo pagados que no
   hacen nada.
2. **RA no puede unirse a un amplicón "plano"**: el adaptador se acopla por
   química click a una modificación propietaria que llevan los primers
   barcodeados de los kits de ONT. Un oligo sintético normal, con la
   secuencia que sea, no la tiene. (Por eso cargar amplicones KYO2 planos
   con RA+ADB produjo cero lecturas: 13 ng = 50 fmol bien cuantificados,
   pero sin adaptador de secuenciación.)

La ruta soportada por ONT para amplicones propios + RA es el **Rapid PCR
Barcoding Kit (SQK-RPB114.24)**, con el enfoque de **cuatro primers**:

```
PCR (una sola reacción por muestra, 4 primers):
  [tag universal F]-[ITS1-F_KYO2]  --->                    (tus primers, baratos)
                                       inserto ITS
                                  <--- [tag universal R]-[ITS2_KYO2]
  +
  primer barcodeado del kit  --> se aposenta sobre el tag universal y añade
  (RLB01..24, propietario)       [extremo de acople para RA]-[barcode ONT]

Después: diluir RA en ADB, acoplar, cargar (idéntico al kit 16S).
```

## Los DOS oligos a pedir

Solo dos — la identidad de muestra ya no va en el primer sino en los
barcodes del kit:

| Nombre | Secuencia 5'→3' | nt |
|---|---|---|
| ITS1-F_KYO2_ONT | `TTTCTGTTGGTGCTGATATTGCTAGAGGAAGTAAAAGTCGTAA` | 43 |
| ITS2_KYO2_ONT | `ACTTGCCTGTCGCTCTATCTTCTTYRCTRCGTTCTTCATC` | 40 |

Estructura: `[tag universal ONT (22 nt)]-[primer KYO2]`. Los tags son los
estándar de ONT para primers con cola:

- Tag F: `TTTCTGTTGGTGCTGATATTGC`
- Tag R: `ACTTGCCTGTCGCTCTATCTTC`

Sin modificaciones químicas, síntesis estándar con desalado basta.
Verifica dímeros/horquillas en OligoAnalyzer antes de pedir, como con
cualquier primer con cola.

## PCR

- Sigue el protocolo del kit RPB114.24 para "custom primers" (enfoque de
  cuatro primers): tus dos primers con cola + el primer barcodeado del kit
  en la misma reacción.
- Los tags no participan en los primeros ciclos: usa el annealing de KYO2
  (~50 °C) los primeros 2–3 ciclos si tu polimerasa lo permite; después el
  amplicón ya lleva las colas y el primer del kit prima sobre ellas.
- Amplicón final ≈ inserto (~260 pb) + primers (39 nt) + tags (44 nt) +
  extensiones barcodeadas del kit ≈ **420–480 pb**.
- Carga: 50 fmol ≈ **13–15 ng** a esa longitud. La cuantificación se hace
  sobre la librería FINAL con RA acoplado.

## Mapa de muestras (barcodes del kit)

| Barcode | Muestra | | Barcode | Muestra |
|---|---|---|---|---|
| RLB01 | S01 | | RLB09 | S09 |
| RLB02 | S02 | | RLB10 | S10 |
| RLB03 | S03 | | RLB11 | S11 |
| RLB04 | S04 | | RLB12 | S12 |
| RLB05 | S05 | | RLB13 | S13 |
| RLB06 | S06 | | RLB14 | S14 |
| RLB07 | S07 | | RLB15 | **MockCom** (Zymo D6300) |
| RLB08 | S08 | | RLB16 | **neg** (blanco) |

Editable en `assets/ont_barcode_map.tsv`. Los Golay dual-index quedan
retirados en este diseño; el pipeline los sigue soportando con
`--demux_mode dual_index` (por defecto) por si se secuencia una librería
del diseño antiguo.

## Análisis

```bash
nextflow run main.nf ... --demux_mode ont_kit --kit SQK-RPB114-24
```

El pipeline entonces demultiplexa con `dorado demux` (barcodes del kit),
recorta tags + primers con cutadapt, y sigue igual (chopper → Emu/UNITE →
phyloseq). La hoja `assets/metadata_template_ont.tsv` da los SampleID de
este modo (S01…S14, MockCom, neg).

---

## Alternativa B: reusar los primers largos Illumina ya sintetizados

Si ya tienes los oligos dual-index largos (P5/P7 + Golay + KYO2), no hace
falta retirarlos: una **segunda PCR de pocos ciclos** injerta los tags
universales sobre los extremos P5/P7 de los amplicones ya barcodeados.

| Nombre | Secuencia 5'→3' | nt |
|---|---|---|
| P5_ONTtag_F | `TTTCTGTTGGTGCTGATATTGCAATGATACGGCGACCACCGAGATCTACAC` | 51 |
| P7_ONTtag_R | `ACTTGCCTGTCGCTCTATCTTCCAAGCAGAAGACGGCATACGAGAT` | 46 |

Estructura: `[tag universal ONT (22 nt)]-[P5 o P7 de Illumina]`. Las
porciones 3' son los primers estándar de amplificación de librerías
Illumina: priman en los extremos exactos del amplicón de la PCR1 y copian
la construcción completa, barcodes Golay incluidos.

Flujo:

1. **PCR1 por muestra** con los primers largos existentes (sin cambios).
2. **Pool** de las 16 muestras + limpieza SPRI 0.8x (los primers/dímeros
   residuales de la PCR1 también se etiquetarían).
3. **PCR2 sobre el pool, 6–10 ciclos**: los dos primers de arriba + UN solo
   primer barcodeado del kit RPB114.24 (p. ej. RLB01) — enfoque de cuatro
   primers sobre el pool. Solo se consume un barcode del kit por corrida;
   la identidad de muestra ya viaja en los Golay. Pocos ciclos: la PCR
   sobre un pool mixto puede generar quimeras entre muestras (aparecerían
   como `invalid_*` en demux_stats.tsv).
4. **Limpieza → Qubit → RA en ADB → carga.** Construcción final ≈ 650–700
   pb; 50 fmol ≈ **21–23 ng**, cuantificados sobre la librería final.

Análisis: `--demux_mode dual_index` (el por defecto) sin ningún cambio —
el demux por Golay ignora la secuencia extra fuera de los barcodes.

## English summary — Alternative B

Reuse the existing long Illumina dual-index primers: PCR1 per sample as
designed, pool + 0.8x SPRI, then a 6–10-cycle PCR2 on the pool with
`TTTCTGTTGGTGCTGATATTGC`+P5 (51 nt) and `ACTTGCCTGTCGCTCTATCTTC`+P7
(46 nt) plus a single RPB114.24 barcoded primer for RA attachment. Load
50 fmol (~21–23 ng at ~650–700 bp). Analyse with the default
`--demux_mode dual_index`; sample identity stays in the Golay indices.

## English summary

The Illumina dual-index constructs are retired. ONT's Rapid Adapter cannot
attach to plain amplicons (the kit's barcoded primers carry a proprietary
click-chemistry modification), so the supported route is the **Rapid PCR
Barcoding Kit (SQK-RPB114.24)** four-primer approach: order only two tailed
primers — `TTTCTGTTGGTGCTGATATTGC` + ITS1-F_KYO2 (43 nt) and
`ACTTGCCTGTCGCTCTATCTTC` + ITS2_KYO2 (40 nt). The kit's RLB01–24 barcoded
primers anneal to those universal tags during the same PCR and add the
barcode plus the RA-attachment end; RA/ADB handling is then identical to
the 16S kit. Final amplicon ≈ 420–480 bp; load 50 fmol (~13–15 ng),
quantified on the FINAL adapted library. Analyse with
`--demux_mode ont_kit --kit SQK-RPB114-24`; sample↔barcode mapping lives in
`assets/ont_barcode_map.tsv` (RLB15 = Zymo D6300 mock, RLB16 = blank).
