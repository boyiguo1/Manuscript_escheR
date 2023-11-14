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

# Session Info ------------------------------------------------------------
# > sessionInfo()
# R version 4.3.0 (2023-04-21)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Ventura 13.6
#
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib
# LAPACK: /Library/Frameworks/R.framework/Versions/4.3-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.11.0
#
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#
# time zone: America/New_York
# tzcode source: internal
#
# attached base packages:
#   [1] stats4    stats     graphics
# [4] grDevices utils     datasets
# [7] methods   base
#
# other attached packages:
#   [1] escheR_1.3.1
# [2] BiocSingular_1.18.0
# [3] ggrepel_0.9.4
# [4] scran_1.30.0
# [5] scater_1.30.0
# [6] scuttle_1.12.0
# [7] TENxPBMCData_1.20.0
# [8] HDF5Array_1.30.0
# [9] rhdf5_2.46.0
# [10] DelayedArray_0.28.0
# [11] SparseArray_1.2.1
# [12] S4Arrays_1.2.0
# [13] abind_1.4-5
# [14] Matrix_1.6-1
# [15] schex_1.16.0
# [16] shiny_1.7.5.1
# [17] ggplot2_3.4.4
# [18] Seurat_5.0.0
# [19] SeuratObject_5.0.0
# [20] sp_2.1-1
# [21] SingleCellExperiment_1.24.0
# [22] SummarizedExperiment_1.32.0
# [23] Biobase_2.62.0
# [24] GenomicRanges_1.54.1
# [25] GenomeInfoDb_1.38.0
# [26] IRanges_2.36.0
# [27] S4Vectors_0.40.1
# [28] BiocGenerics_0.48.1
# [29] MatrixGenerics_1.14.0
# [30] matrixStats_1.1.0
# [31] igraph_1.5.1
#
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.0-3
# [2] bitops_1.0-7
# [3] httr_1.4.7
# [4] RColorBrewer_1.1-3
# [5] doParallel_1.0.17
# [6] tools_4.3.0
# [7] sctransform_0.4.1
# [8] backports_1.4.1
# [9] utf8_1.2.4
# [10] R6_2.5.1
# [11] DT_0.30
# [12] uwot_0.1.16
# [13] lazyeval_0.2.2
# [14] rhdf5filters_1.14.1
# [15] withr_2.5.2
# [16] gridExtra_2.3
# [17] progressr_0.14.0
# [18] cli_3.6.1
# [19] spatstat.explore_3.2-5
# [20] fastDummies_1.7.3
# [21] entropy_1.3.1
# [22] sass_0.4.7
# [23] spatstat.data_3.0-3
# [24] ggridges_0.5.4
# [25] pbapply_1.7-2
# [26] Rsamtools_2.18.0
# [27] parallelly_1.36.0
# [28] sessioninfo_1.2.2
# [29] attempt_0.3.1
# [30] maps_3.4.1.1
# [31] limma_3.58.1
# [32] rstudioapi_0.15.0
# [33] RSQLite_2.3.3
# [34] generics_0.1.3
# [35] BiocIO_1.12.0
# [36] spatstat.random_3.2-1
# [37] ica_1.0-3
# [38] car_3.1-2
# [39] dplyr_1.1.3
# [40] ggbeeswarm_0.7.2
# [41] fansi_1.0.5
# [42] lifecycle_1.0.4
# [43] yaml_2.3.7
# [44] edgeR_4.0.1
# [45] carData_3.0-5
# [46] BiocFileCache_2.10.1
# [47] Rtsne_0.16
# [48] paletteer_1.5.0
# [49] grid_4.3.0
# [50] blob_1.2.4
# [51] dqrng_0.3.1
# [52] promises_1.2.1
# [53] ExperimentHub_2.10.0
# [54] crayon_1.5.2
# [55] miniUI_0.1.1.1
# [56] lattice_0.22-5
# [57] beachmat_2.18.0
# [58] cowplot_1.1.1
# [59] KEGGREST_1.42.0
# [60] magick_2.8.1
# [61] metapod_1.10.0
# [62] pillar_1.9.0
# [63] rjson_0.2.21
# [64] future.apply_1.11.0
# [65] codetools_0.2-19
# [66] leiden_0.4.3
# [67] glue_1.6.2
# [68] data.table_1.14.8
# [69] vctrs_0.6.4
# [70] png_0.1-8
# [71] spam_2.10-0
# [72] gtable_0.3.4
# [73] rematch2_2.1.2
# [74] cachem_1.0.8
# [75] mime_0.12
# [76] survival_3.5-7
# [77] iterators_1.0.14
# [78] fields_15.2
# [79] bluster_1.12.0
# [80] statmod_1.5.0
# [81] interactiveDisplayBase_1.40.0
# [82] ellipsis_0.3.2
# [83] fitdistrplus_1.1-11
# [84] ROCR_1.0-11
# [85] nlme_3.1-163
# [86] bit64_4.0.5
# [87] filelock_1.0.2
# [88] RcppAnnoy_0.0.21
# [89] bslib_0.5.1
# [90] irlba_2.3.5.1
# [91] vipor_0.4.5
# [92] KernSmooth_2.23-22
# [93] colorspace_2.1-0
# [94] DBI_1.1.3
# [95] tidyselect_1.2.0
# [96] bit_4.0.5
# [97] compiler_4.3.0
# [98] curl_5.1.0
# [99] BiocNeighbors_1.20.0
# [100] plotly_4.10.3
# [101] rtracklayer_1.62.0
# [102] scales_1.2.1
# [103] hexbin_1.28.3
# [104] lmtest_0.9-40
# [105] rappdirs_0.3.3
# [106] goftest_1.2-3
# [107] stringr_1.5.0
# [108] SpatialExperiment_1.12.0
# [109] digest_0.6.33
# [110] spatstat.utils_3.0-4
# [111] benchmarkmeData_1.0.4
# [112] XVector_0.42.0
# [113] htmltools_0.5.7
# [114] pkgconfig_2.0.3
# [115] sparseMatrixStats_1.14.0
# [116] dbplyr_2.4.0
# [117] fastmap_1.1.1
# [118] rlang_1.1.2
# [119] htmlwidgets_1.6.2
# [120] spatialLIBD_1.13.4
# [121] DelayedMatrixStats_1.24.0
# [122] farver_2.1.1
# [123] jquerylib_0.1.4
# [124] zoo_1.8-12
# [125] jsonlite_1.8.7
# [126] BiocParallel_1.36.0
# [127] config_0.3.2
# [128] RCurl_1.98-1.13
# [129] magrittr_2.0.3
# [130] GenomeInfoDbData_1.2.11
# [131] dotCall64_1.1-0
# [132] patchwork_1.1.3
# [133] Rhdf5lib_1.24.0
# [134] munsell_0.5.0
# [135] Rcpp_1.0.11
# [136] viridis_0.6.4
# [137] reticulate_1.34.0
# [138] stringi_1.7.12
# [139] zlibbioc_1.48.0
# [140] MASS_7.3-60
# [141] AnnotationHub_3.10.0
# [142] plyr_1.8.9
# [143] parallel_4.3.0
# [144] listenv_0.9.0
# [145] deldir_1.0-9
# [146] Biostrings_2.70.1
# [147] splines_4.3.0
# [148] tensor_1.5
# [149] locfit_1.5-9.8
# [150] ggpubr_0.6.0
# [151] spatstat.geom_3.2-7
# [152] ggsignif_0.6.4
# [153] RcppHNSW_0.5.0
# [154] reshape2_1.4.4
# [155] ScaledMatrix_1.10.0
# [156] BiocVersion_3.18.0
# [157] XML_3.99-0.15
# [158] golem_0.4.1
# [159] BiocManager_1.30.22
# [160] tweenr_2.0.2
# [161] foreach_1.5.2
# [162] httpuv_1.6.12
# [163] polyclip_1.10-6
# [164] RANN_2.6.1
# [165] tidyr_1.3.0
# [166] purrr_1.0.2
# [167] future_1.33.0
# [168] benchmarkme_1.0.8
# [169] scattermore_1.2
# [170] ggforce_0.4.1
# [171] rsvd_1.0.5
# [172] broom_1.0.5
# [173] xtable_1.8-4
# [174] restfulr_0.0.15
# [175] RSpectra_0.16-1
# [176] rstatix_0.7.2
# [177] later_1.3.1
# [178] viridisLite_0.4.2
# [179] tibble_3.2.1
# [180] memoise_2.0.1
# [181] beeswarm_0.4.0
# [182] AnnotationDbi_1.64.1
# [183] GenomicAlignments_1.38.0
# [184] cluster_2.1.4
# [185] shinyWidgets_0.8.0
# [186] globals_0.16.2
# [187] concaveman_1.1.0
