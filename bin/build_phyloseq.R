#!/usr/bin/env Rscript

# Construye un objeto phyloseq desde la salida combinada de Emu, con QC
# básico.
#
#   build_phyloseq.R <combined-species-counts.tsv> <metadata.tsv> [outdir] [min_reads]
#
# Los metadatos necesitan una columna SampleID que coincida con los nombres
# de muestra de la hoja de barcodes (F1R1, F1R2, ... en assets/barcodes.tsv),
# que son los nombres de columna en la tabla combinada de Emu.

suppressPackageStartupMessages({
  library(phyloseq)
  library(readr)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("uso: build_phyloseq.R <counts.tsv> <metadata.tsv> [outdir] [min_reads]")
}

counts_file <- args[1]
meta_file   <- args[2]
outdir      <- if (length(args) >= 3) args[3] else "phyloseq_out"
min_reads   <- if (length(args) >= 4) as.numeric(args[4]) else 1000

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# UNITE usa kingdom (Fungi) donde GTDB/SILVA usan superkingdom; se aceptan
# ambos. "clade" aparece en algunas construcciones de UNITE.
RANKS <- c("kingdom", "superkingdom", "clade", "phylum", "class", "order",
           "family", "genus", "species")

# ---------------------------------------------------------------- entradas

raw <- read_tsv(counts_file, show_col_types = FALSE)
meta <- read_tsv(meta_file, show_col_types = FALSE)

if (!"SampleID" %in% names(meta)) {
  stop("los metadatos deben contener una columna SampleID")
}

# Las columnas de muestra son todo lo que no sea taxonomía: los nombres
# vienen de la hoja de barcodes, así que no hay un prefijo fijo que buscar.
non_sample <- c("tax_id", "lineage", RANKS)
sample_cols <- setdiff(names(raw), non_sample)
if (!length(sample_cols)) {
  stop("no se encontraron columnas de muestra en ", counts_file)
}

missing <- setdiff(sample_cols, meta$SampleID)
if (length(missing)) {
  message("ADVERTENCIA: sin fila de metadatos para: ", paste(missing, collapse = ", "))
}

# ---------------------------------------------------------------- construcción

# Los tax_id de Emu pueden repetirse entre linajes distintos en algunas
# construcciones, así que usa IDs posicionales de fila en lugar de tax_id
# como clave de taxones.
ids <- paste0("t", seq_len(nrow(raw)))

otu <- as.matrix(raw[, sample_cols])
otu[is.na(otu)] <- 0
# El EM de Emu produce estimaciones fraccionarias (p. ej. 39879.125).
# DESeq2 y rarefy_even_depth() requieren enteros.
otu <- round(otu)
rownames(otu) <- ids

tax <- as.matrix(raw[, intersect(RANKS, names(raw))])
tax[tax == "" | tax == "-"] <- NA
rownames(tax) <- ids

md <- tibble(SampleID = sample_cols) |>
  left_join(meta, by = "SampleID") |>
  as.data.frame()
rownames(md) <- md$SampleID

ps_all <- phyloseq(
  otu_table(otu, taxa_are_rows = TRUE),
  tax_table(tax),
  sample_data(md)
)

# ------------------------------------------------------- QC de profundidad

depth <- tibble(SampleID = sample_names(ps_all), reads = sample_sums(ps_all)) |>
  arrange(desc(reads)) |>
  mutate(retained = reads >= min_reads)

write_tsv(depth, file.path(outdir, "sample_depth.tsv"))

# Las muestras de baja profundidad producen abundancias relativas
# inestables Y riqueza aparente inflada. Excluirlas no es opcional si
# planeas ordenar. El control negativo debe caer aquí; si no cae,
# investiga contaminación antes de interpretar nada.
dropped <- depth |> filter(!retained) |> pull(SampleID)
if (length(dropped)) {
  message("Excluyendo ", length(dropped), " muestra(s) por debajo de ", min_reads,
          " lecturas: ", paste(dropped, collapse = ", "))
}

ps <- prune_samples(!(sample_names(ps_all) %in% dropped), ps_all)
ps <- prune_taxa(taxa_sums(ps) > 0, ps)

saveRDS(ps, file.path(outdir, "phyloseq.rds"))
saveRDS(ps_all, file.path(outdir, "phyloseq_unfiltered.rds"))

capture.output(print(ps), file = file.path(outdir, "phyloseq_summary.txt"))

# ---------------------------------------------------------------- tablas

ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

# UNITE etiqueta los taxones sin resolver con placeholders tipo
# "unidentified" o hipótesis de especie "Genus_sp". En un sistema poco
# estudiado esta fracción es un hallazgo, no un error.
sp <- as.character(tax_table(ps)[, "species"])
placeholder <- grepl("unidentified|_sp$|_sp\\.|Incertae", sp, ignore.case = TRUE) & !is.na(sp)

if (any(placeholder)) {
  tibble(
    SampleID        = sample_names(ps_rel),
    uncharacterised = sample_sums(prune_taxa(placeholder, ps_rel))
  ) |>
    arrange(desc(uncharacterised)) |>
    write_tsv(file.path(outdir, "uncharacterised_fraction.tsv"))
}

# Shannon y Simpson son mucho más confiables que la riqueza Observada sobre
# salida de EM, que distribuye las lecturas entre referencias similares.
alpha <- estimate_richness(ps, measures = c("Shannon", "Simpson", "InvSimpson")) |>
  rownames_to_column("SampleID")
write_tsv(alpha, file.path(outdir, "alpha_diversity.tsv"))

# ---------------------------------------------------------------- gráficos

ok <- try({
  ps_gen <- tax_glom(ps_rel, "genus", NArm = FALSE)
  n_top  <- min(15, ntaxa(ps_gen))
  top    <- names(sort(taxa_sums(ps_gen), decreasing = TRUE))[seq_len(n_top)]

  p <- plot_bar(prune_taxa(top, ps_gen), fill = "genus") +
    labs(title = paste("Top", n_top, "géneros"), y = "Abundancia relativa") +
    theme_bw() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 90, hjust = 1, size = 6))
  ggsave(file.path(outdir, "genus_barplot.pdf"), p, width = 12, height = 8)
}, silent = TRUE)
if (inherits(ok, "try-error")) message("gráfico de barras por género omitido: ", ok)

# Emu no produce filogenia, así que UniFrac no está disponible. Se usa
# Bray-Curtis.
ok <- try({
  if (nsamples(ps_rel) >= 4) {
    ord <- ordinate(ps_rel, "PCoA", "bray")
    p <- plot_ordination(ps_rel, ord) +
      geom_text(aes(label = sample_names(ps_rel)), size = 2.5, vjust = -1) +
      labs(title = "PCoA (Bray-Curtis)") + theme_bw()
    ggsave(file.path(outdir, "ordination_pcoa.pdf"), p, width = 9, height = 7)
  }
}, silent = TRUE)
if (inherits(ok, "try-error")) message("ordenación omitida: ", ok)

# ---------------------------------------------------------------- estadística

ok <- try({
  if (requireNamespace("vegan", quietly = TRUE) && nsamples(ps_rel) >= 6) {
    d  <- phyloseq::distance(ps_rel, "bray")
    smd <- data.frame(sample_data(ps_rel))
    vars <- setdiff(names(smd), "SampleID")
    vars <- vars[vapply(smd[vars], function(x) {
      k <- length(unique(x[!is.na(x)]))
      k > 1 && k < nrow(smd)
    }, logical(1))]

    for (v in vars) {
      res <- try(vegan::adonis2(as.formula(paste("d ~", v)),
                                data = smd, permutations = 999), silent = TRUE)
      if (!inherits(res, "try-error")) {
        capture.output(print(res),
                       file = file.path(outdir, paste0("adonis_", v, ".txt")))
      }
    }
  }
}, silent = TRUE)

capture.output(sessionInfo(), file = file.path(outdir, "sessionInfo.txt"))
message("salida de phyloseq escrita en: ", outdir)
