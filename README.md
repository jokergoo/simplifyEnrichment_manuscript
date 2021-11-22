# simplifyEnrichment_manuscript

Scripts for the analysis in [simplifyEnrichment](https://github.com/jokergoo/simplifyEnrichment) manuscript:

- `collect_sim_mat.R`: It performes gene set enrichment analysis on [EBI Expression Atlas datasets](https://www.ebi.ac.uk/gxa/download) with various ontologies. It also calculates similarity matrices with different measurements on the significant terms.
- `run_example_random_GO.R`: It runs various clustering methods on the similarity matrices that were generated from random GO terms.
- `run_examples_EBI.R`: It runs various clustering methods on the similarity matrices from [EBI Expression Atlas datasets](https://www.ebi.ac.uk/gxa/download).
- `compare_semantic_and_overlap.R`: It compares clusterings from different similarity measurements, i.e., semantic similarity and gene overlap similarities.
- `test_partition_methods.R`: It tests the effect of partitionning methods in the binary cut clustering with rando GO datasets.
- `website.R`: It deploys all the results to GitHub Page.
- `figures/`: Script for making figures in the manuscript.

Many analysis in the scripts were sent to the computing cluster with [the bsub package](https://github.com/jokergoo/bsub).

### Session info for the analysis

```
R version 4.1.0 (2021-05-18)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: CentOS Linux 7 (Core)

Matrix products: default
BLAS/LAPACK: /usr/lib64/libopenblas-r0.3.3.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C
 [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8
 [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8
 [7] LC_PAPER=en_US.UTF-8       LC_NAME=C
 [9] LC_ADDRESS=C               LC_TELEPHONE=C
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C

attached base packages:
 [1] stats4    grid      parallel  stats     graphics  grDevices utils
 [8] datasets  methods   base

other attached packages:
 [1] gt_0.3.1                 BiocManager_1.30.16      flexclust_1.4-0
 [4] modeltools_0.2-23        lattice_0.20-44          reactome.db_1.76.0
 [7] DO.db_2.9                AnnotationDbi_1.54.1     IRanges_2.26.0
[10] S4Vectors_0.30.0         Biobase_2.52.0           gridGraphics_0.5-1
[13] testthat_3.0.4           gridExtra_2.3            apcluster_1.4.8
[16] mclust_5.4.7             cowplot_1.1.1            ggplot2_3.3.5
[19] knitr_1.34               cola_1.99.4              msigdbr_7.4.1
[22] DOSE_3.18.2              dynamicTreeCut_1.63-1    dbscan_1.1-8
[25] MCL_1.0                  igraph_1.2.6             clusterProfiler_4.2.0
[28] simplifyEnrichment_1.5.1 BiocGenerics_0.38.0      colorout_1.2-2

loaded via a namespace (and not attached):
  [1] shadowtext_0.0.9       circlize_0.4.13        fastmatch_1.1-3
  [4] plyr_1.8.6             lazyeval_0.2.2         proxyC_0.2.1
  [7] splines_4.1.0          BiocParallel_1.26.2    GenomeInfoDb_1.28.4
 [10] digest_0.6.27          htmltools_0.5.2        foreach_1.5.1
 [13] yulab.utils_0.0.2      GOSemSim_2.18.1        viridis_0.6.1
 [16] GO.db_3.13.0           fansi_0.5.0            magrittr_2.0.1
 [19] memoise_2.0.0          tm_0.7-8               cluster_2.1.2
 [22] doParallel_1.0.16      ComplexHeatmap_2.8.0   Biostrings_2.60.2
 [25] graphlayouts_0.7.1     RcppParallel_5.1.4     matrixStats_0.61.0
 [28] enrichplot_1.12.2      colorspace_2.0-2       blob_1.2.2
 [31] ggrepel_0.9.1          xfun_0.24              dplyr_1.0.7
 [34] microbenchmark_1.4-7   crayon_1.4.1           RCurl_1.98-1.5
 [37] jsonlite_1.7.2         scatterpie_0.1.7       impute_1.66.0
 [40] brew_1.0-6             iterators_1.0.13       ape_5.5
 [43] glue_1.4.2             polyclip_1.10-0        gtable_0.3.0
 [46] zlibbioc_1.38.0        XVector_0.32.0         GetoptLong_1.0.5
 [49] shape_1.4.6            scales_1.1.1           DBI_1.1.1
 [52] Rcpp_1.0.7             viridisLite_0.4.0      clue_0.3-59
 [55] tidytree_0.3.5         bit_4.0.4              httr_1.4.2
 [58] fgsea_1.18.0           RColorBrewer_1.1-2     ellipsis_0.3.2
 [61] pkgconfig_2.0.3        farver_2.1.0           utf8_1.2.1
 [64] ggplotify_0.1.0        tidyselect_1.1.1       rlang_0.4.12
 [67] reshape2_1.4.4         munsell_0.5.0          tools_4.1.0
 [70] cachem_1.0.5           downloader_0.4         generics_0.1.0
 [73] RSQLite_2.2.8          stringr_1.4.0          fastmap_1.1.0
 [76] ggtree_3.0.4           org.Hs.eg.db_3.13.0    babelgene_21.4
 [79] bit64_4.0.5            tidygraph_1.2.0        purrr_0.3.4
 [82] KEGGREST_1.32.0        ggraph_2.0.5           nlme_3.1-152
 [85] slam_0.1-48            aplot_0.1.0            xml2_1.3.2
 [88] compiler_4.1.0         rstudioapi_0.13        png_0.1-7
 [91] treeio_1.16.2          tibble_3.1.4           tweenr_1.0.2
 [94] stringi_1.7.4          Matrix_1.3-4           markdown_1.1
 [97] vctrs_0.3.8            pillar_1.6.2           lifecycle_1.0.0
[100] eulerr_6.1.1           GlobalOptions_0.1.2    irlba_2.3.3
[103] data.table_1.14.0      bitops_1.0-7           patchwork_1.1.1
[106] qvalue_2.24.0          R6_2.5.1               codetools_0.2-18
[109] MASS_7.3-54            assertthat_0.2.1       rjson_0.2.20
[112] withr_2.4.2            GenomeInfoDbData_1.2.6 expm_0.999-6
[115] ggfun_0.0.4            class_7.3-19           tidyr_1.1.3
[118] skmeans_0.2-13         Cairo_1.5-12.2         ggforce_0.3.3
[121] NLP_0.2-1
```

