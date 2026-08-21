# Changelog

## 1.1.0 — 2026-08-21

- **Nuevo modo de demultiplexado `ont_kit`** para el rediseño compatible
  con Rapid Adapter: primers KYO2 con tags universales ONT
  (`TTTCTGTTGGTGCTGATATTGC` / `ACTTGCCTGTCGCTCTATCTTC`) + barcodes del
  Rapid PCR Barcoding Kit (SQK-RPB114.24). `dorado demux` + recorte de
  tags/primers con cutadapt (`--revcomp --discard-untrimmed`). Ver
  `docs/PRIMERS_ONT_RA.md` para el rediseño y su justificación (RA no se
  acopla a amplicones sin la modificación propietaria de los primers del
  kit — los amplicones planos cargados con RA/ADB producen cero lecturas).
- Nuevos assets: `ont_barcode_map.tsv` (RLB01–16; 15 = mock Zymo D6300,
  16 = blanco) y `metadata_template_ont.tsv`.
- El modo `dual_index` original queda como valor por defecto, sin cambios.

## 1.0.0 — 2026-08-19

Primera versión, derivada de ont16s 1.0.0. Cambios respecto a ont16s:

- **Demultiplexado**: `dorado demux` (kit SQK-16S114-24) reemplazado por
  `bin/demux_dual_barcodes.py` — dos rondas de cutadapt sobre la hoja
  `assets/barcodes.tsv` (primers ITS1-F_KYO2 / ITS2_KYO2, 4 barcodes Golay
  F × 4 R = 16 combinaciones; sample 15 = MockCom, 16 = control negativo).
  Ambos barcodes obligatorios; combinaciones fuera de la hoja contadas como
  index hopping. Validado con lecturas sintéticas al 3% de error: 95.8%
  asignadas, 0.2% mal asignadas, 0/20 lecturas basura aceptadas.
- **Filtrado**: ventana 150–600 pb (amplicón KYO2 ~300 pb, con margen para
  la variación de ITS1) en lugar de 1200–1800; `min_reads` 1000 en lugar
  de 5000.
- **Taxonomía**: bases preconstruidas UNITE de Emu (`unite-fungi`,
  `unite-all`) en lugar de GTDB/SILVA; `setup_databases.sh` reescrito.
- **phyloseq**: detección de columnas de muestra por la hoja de barcodes en
  lugar del prefijo `barcode*`; rangos con `kingdom` (UNITE); placeholders
  "unidentified"/SH en lugar de los de GTDB.
- **Entornos**: cutadapt añadido al entorno `ont`.
