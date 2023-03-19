# Packages ----------------------------------------------------------------
library(igraph)
library(schex)
library(TENxPBMCData)
library(scater)
library(scran)
library(ggrepel)



# Code --------------------------------------------------------------------
# This section of code was directly taken from schex vignettes.
# https://github.com/SaskiaFreytag/schex/blob/master/vignettes/using_schex.Rmd

# Data importing and preprocessing
tenx_pbmc3k <- TENxPBMCData(dataset = "pbmc3k")

rownames(tenx_pbmc3k) <- uniquifyFeatureNames(rowData(tenx_pbmc3k)$ENSEMBL_ID,
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
tenx_pbmc3k <- runPCA(tenx_pbmc3k)

# Changed the seed number from 10 to 100
set.seed(100)

tenx_pbmc3k <- runUMAP(tenx_pbmc3k, dimred = "PCA", spread = 1,
                       min_dist = 0.4)
snn_gr <- buildSNNGraph(tenx_pbmc3k, use.dimred = "PCA", k = 50)
clusters <- cluster_louvain(snn_gr)
tenx_pbmc3k$cluster <- factor(clusters$membership)

tenx_pbmc3k <- make_hexbin(tenx_pbmc3k, nbins = 40,
                           dimension_reduction = "UMAP", use_dims=c(1,2))

plot_hexbin_density(tenx_pbmc3k)

plot_hexbin_meta(tenx_pbmc3k, col="cluster", action="majority")
plot_hexbin_meta(tenx_pbmc3k, col="total", action="median")


gene_id <-"POMGNT1"

orginin_UMAP_P <- plot_hexbin_feature_plus(
  tenx_pbmc3k,
  col="cluster", type="logcounts",
  feature="POMGNT1", action="mean") +
  labs(title = "") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_blank())

ggsave("figs/schex_UMAP.pdf",
       orginin_UMAP_P,
       height = 7,
       width = 7)

# END (Code from schex vignettes ends here)

# Multi-dimensional UMAP --------------------------------------------------
# New code created by Boyi Guo to show multi-dimensional UMAP
improved_UMP_P <- plot_hexbin_feature_plus(
  tenx_pbmc3k,
  col="cluster", type="logcounts",
  feature="POMGNT1", action="mean") +
  geom_hex(aes(color = meta),
           stat = "identity",
           fill = "transparent", size =0.5)

# Remove convex hulls
improved_UMP_P$layers <- improved_UMP_P$layers[-2]

improved_UMP_P <- improved_UMP_P +
  labs(title = "") +
  theme_classic(base_size = 20) +
  theme(legend.title = element_blank())

ggsave("figs/escher_UMAP.pdf",improved_UMP_P,
       height = 7,
       width = 7)

