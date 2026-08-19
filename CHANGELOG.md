# Changelog

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
