# Packages ----------------------------------------------------------------
library(SpatialExperiment)
library(spatialLIBD)
library(here)
library(tidyverse)
library(escheR)             # Proposed package



# Helper Function ---------------------------------------------------------

# Function to calculate the gene_name1, gene_name2 co-localization
# of an SPE object
# Return a SPE object with co-localization status as one variable in the colData
calc_coloc <- function(spe, gene_name1, gene_name2, sample_id = NULL) {
  if (!is.null(sample_id)) {
    if (!sample_id %in% unique(spe$sample_id)) {
      stop("Can't find sample_id")
    }
    spe <- spe[, spe$sample_id == sample_id]
  }

  if (!all(c(gene_name1, gene_name2) %in% rowData(spe)$gene_name)) {
    stop("Can't find gene names")
  }

  gene1_idx <- which(rowData(spe)$gene_name == gene_name1)
  gene2_idx <- which(rowData(spe)$gene_name == gene_name2)

  cnt_df <- counts(spe)[c(gene1_idx, gene2_idx), ] |>
    as.matrix() |>
    t() |>
    data.frame()

  stopifnot(setequal(names(cnt_df), rowData(spe)[names(cnt_df), "gene_id"]))

  old_names <- names(cnt_df)
  names(cnt_df) <- paste0("gene", 1:2)
  # Check if the correct rows are selected

  coloc_cat <- cnt_df |>
    pmap_chr(.f = function(gene1, gene2) {
      if (gene1 != 0 && gene2 != 0) {
        return(paste0(gene_name1, "&", gene_name2))
      } else if (gene1 != 0) {
        return(gene_name1)
      } else if (gene2 != 0) {
        return(gene_name2)
      } else {
        return("Neither")
      }
    })

  colData(spe)$coloc <- coloc_cat |> factor()

  return(spe)
}


# Code --------------------------------------------------------------------

if (!exists("spe"))
  spe <- fetch_data("spatialDLPFC_Visium")

# Subset sample
spe <- spe[, spe$sample_id == "Br8667_mid"]


# Create colocalization variable in colData(spe)
spe <- calc_coloc(spe,
                  gene_name1 = "EFNA5",
                  gene_name2 = "EPHA5")

# Rename spatial domain
colData(spe)$spDomain <- factor(
  spe$BayesSpace_harmony_09,
  levels = c(1, 2, 3, 5, 8, 4, 7, 6, 9),
  labels = c("L1", "L1", "L2", "L3", "L4", "L5", "L6", "WM", "WM")
)


# Colocalization plot (Panel A) -------------------------------------------
p_coloc <- spatialLIBD::vis_clus(
  spe = spe,
  clustervar = "coloc",
  spatial = FALSE,
  # ... = " LIBD Layers",
  # point_size = 10
) +
  scale_fill_manual(
    name = "",
    limits = c("EFNA5&EPHA5",
               "EFNA5",
               "EPHA5",
               "Neither"),
    values = c("EFNA5&EPHA5" = "#483838",
               "EFNA5" = "#96E5D1",
               "EPHA5" = "#89CFFD",
               "Neither" = "#CCCCCC40"
    )
  ) +
  labs(title = "") +
  theme(legend.position = "none")


ggsave(
  here("figs", "coloc.pdf"),
  plot = p_coloc,
  height = 7,
  width = 7)

# Spatial Domain Plot (Panel B) --------------------------------------------
p_layer <- spatialLIBD::vis_clus(
  spe = spe,
  clustervar = "spDomain",
  spatial = FALSE,
  colors = libd_layer_colors
  # ... = " LIBD Layers",
  # point_size = 10
) +
  labs(title = "")+
  theme(legend.position = "none")

ggsave(
  here("figs", "spd_plot.pdf"),
  plot = p_layer,
  height = 7,
  width = 7
)

# Spatial Domain Plot with Escher (Panel C) --------------------------------
p_layer_esher <- (
  makeCOCOON(spe) |>
    add_shade.SpatialExperiment(
      var = "spDomain")
) +
  scale_color_manual(
    name = "",
    values = libd_layer_colors |>
      setNames(c(paste0("L", 1:6), "WM", "NA", "WM2"))
  )  +
  labs(title = "") +
  theme(legend.position = "none")

ggsave(
  here("figs, spd_escheR.pdf"),
  plot = p_layer_esher,
  height = 7,
  width = 7)

# Multi-dimensional Plot (Panel D) --------------------------------------
p_coloc_esher <- (p_coloc |>
                    add_shade.SpatialExperiment(var = "spDomain")) +
  scale_color_manual(
    name = "",
    values = libd_layer_colors |>
      setNames(c(paste0("L", 1:6), "WM", "NA", "WM2"))
  )  +
  labs(title = "") +
  theme(legend.position = "none")

ggsave(
  here("figs", "coloc_escher.pdf"),
  plot = p_coloc_esher,
  height = 7,
  width = 7)

