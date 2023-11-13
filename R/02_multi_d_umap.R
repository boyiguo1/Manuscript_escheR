# Packages ----------------------------------------------------------------
library(igraph)
library(schex)
library(TENxPBMCData)
library(scater)
library(scran)
library(scuttle)
library(ggrepel)
library(BiocSingular)
library(escheR)
stopifnot(package.version("escheR") >= "1.3.1")

if(Biobase::package.version("Matrix") > "1.6-1")
  # Bug report on Matrix 1.6-2: https://github.com/bwlewis/irlba/issues/70
  stop("Matrix 1.6-2 may not function properly to calculate PCA.")



# Data Reprocessing Code --------------------------------------------------------------------
# This section of code was directly taken from schex vignettes.
# https://github.com/SaskiaFreytag/schex/blob/master/vignettes/using_schex.Rmd
set.seed(3)
# Data importing and preprocessing
tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")

rownames(tenx_pbmc3k) <- scuttle::uniquifyFeatureNames(rowData(tenx_pbmc3k)$ENSEMBL_ID,
                                                       rowData(tenx_pbmc3k)$Symbol_TENx)
rowData(tenx_pbmc3k)$Mito <- grepl("^MT-", rownames(tenx_pbmc3k))
colData(tenx_pbmc3k) <- cbind(colData(tenx_pbmc3k),
                              perCellQCMetrics(tenx_pbmc3k,
                                               subsets=list(Mt=rowData(tenx_pbmc3k)$Mito)))
rowData(tenx_pbmc3k) <- cbind(rowData(tenx_pbmc3k),
                              perFeatureQCMetrics(tenx_pbmc3k))

tenx_pbmc3k <- tenx_pbmc3k[, !colData(tenx_pbmc3k)$subsets_Mt_percent > 50]

libsize_drop <- isOutlier(tenx_pbmc3k$total,
                          nmads = 3,type = "lower", log = TRUE)
feature_drop <- isOutlier(tenx_pbmc3k$detected,
                          nmads = 3, type = "lower", log = TRUE)

tenx_pbmc3k <- tenx_pbmc3k[, !(libsize_drop | feature_drop)]
rm_ind <- calculateAverage(tenx_pbmc3k)<0
tenx_pbmc3k <- tenx_pbmc3k[!rm_ind,]
tenx_pbmc3k <- scater::logNormCounts(tenx_pbmc3k)
tenx_pbmc3k <- scater::runPCA(tenx_pbmc3k)

# Changed the seed number from 10 to 100

#
tenx_pbmc3k <- runUMAP(tenx_pbmc3k, dimred = "PCA", spread = 1,
                       min_dist = 0.4)
snn_gr <- buildSNNGraph(tenx_pbmc3k, use.dimred = "PCA", k = 50)
clusters <- cluster_louvain(snn_gr)
tenx_pbmc3k$cluster <- factor(clusters$membership)


# schex Implementation ----------------------------------------------------
set.seed(1)
tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, nbins = 50,
                           dimension_reduction = "UMAP", use_dims=c(1,2))

# plot_hexbin_density(tenx_pbmc3k)

# plot_hexbin_meta(tenx_pbmc3k, col="cluster", action="majority")
# plot_hexbin_meta(tenx_pbmc3k, col="total", action="median")

orginin_UMAP_P <- plot_hexbin_feature_plus(
  tenx_pbmc3k,
  col="cluster", type="logcounts",
  feature="POMGNT1", action="mean") +
  geom_segment(
    x = 0.5, xend = -3,
    y = -1, yend = -1,
    linewidth = 2,
    arrow = arrow(length = unit(0.3, "inches"))
  ) +
  geom_segment(
    x = 2.5, xend = 0,
    y = 2.5, yend = 4,
    linewidth = 2,
    arrow = arrow(length = unit(0.3, "inches"))
  ) +
  scale_fill_viridis_c(
    name = "POMGNT1",
    limits = c(0,1.2)
    ) +
  labs(
    title = expression(italic("schex")),
    x = "UMAP 1",
    y = "UMAP 2")  +
  theme_classic(base_size = 20) +
 theme(legend.position="none",
       plot.title = element_text(hjust = 0.5))
# END (Code from schex vignettes ends here)

# escheR Visualization--------------------------------------------------
# New code created by Boyi Guo to show multi-dimensional UMAP
tenx_pbmc3k$logcounts_POMGNT1 <-logcounts(tenx_pbmc3k)["POMGNT1",]
set.seed(1)
escheR_plot <- (make_escheR(tenx_pbmc3k, dimred = "UMAP") |>
                  add_ground_bin("cluster", bins = 50) |>
                  add_fill_bin("logcounts_POMGNT1", bins = 50, fun = "mean")) +
  geom_segment(
    x = 0.5, xend = -3,
    y = -1, yend = -1,
    linewidth = 2,
    arrow = arrow(length = unit(0.3, "inches"))
  ) +
  geom_segment(
    x = 2.5, xend = 0,
    y = 2.5, yend = 4,
    linewidth = 2,
    arrow = arrow(length = unit(0.3, "inches"))
  ) +
  scale_fill_viridis_c(
    name = "POMGNT1",
    limits = c(0,1.2)
    ) +
  scale_color_discrete(name = "Cluster") +
  labs(
    title = expression(italic("escheR")),
    x = "UMAP 1",
    y = "UMAP 2"
  )  +
  theme_classic(base_size = 20) +
  theme(legend.text=element_text(size=12),
        plot.title = element_text(hjust = 0.5))
#


# Final Plot --------------------------------------------------------------
p_fig2 <- ggpubr::ggarrange(
  orginin_UMAP_P,
  escheR_plot,
  nrow = 1,
  labels = "AUTO",
  font.label = list(size = 20, color = "black", face = "bold"),
  legend = "bottom",
  # common.legend = TRUE,
  legend.grob = ggpubr::get_legend(escheR_plot, position = "bottom")
)


ggsave("Manuscript/figure/embedding.tiff", p_fig2,
       height = 8,
       width = 12)

