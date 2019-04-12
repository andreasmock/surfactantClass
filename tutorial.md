Predicting the surfactant class from 450k arrays
================
Andreas Mock and Kolja Pocher

#### Introduction

The following tutorial describes the prediction of the"surfactant high" and "surfactant low" classes from methylome data (450k arrays). The predictor requires the CellMix R package.

``` r
library(CellMix)
```

The predictor can be loaded by:

``` r
load("https://github.com/andreasmock/surfactantClass/surfactantClass.RData")
```

#### Download example data

To illustrate the functionality of the predictor, we download a examplary lung adenocarcinoma sample analysed with the Illumina 450k array from an independent publish study using the Gene Expression Omnibus (Bjaanæs et al., 2016, [GSE66836](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?&acc=GSE66836)). The sample id is GSM1632880.

``` r
sample <- read.table(file="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?mode=raw&is_datatable=true&acc=GSM1632880&id=21811&db=GeoDb_blob124", sep="\t", header = TRUE)
```

The object `sample` comprises the beta values of 456946 CpGs.

``` r
dim(sample)
```

    ## [1] 456946      2

``` r
head(sample)
```

    ##       ID_REF     VALUE
    ## 1 cg00000029 0.2034067
    ## 2 cg00000108 0.8475302
    ## 3 cg00000109 0.8869764
    ## 4 cg00000165 0.1566855
    ## 5 cg00000236 0.8175950
    ## 6 cg00000289 0.1943085

#### Predict surfactant class

The function `surfactantClass` applies a deconvolution algorithm to estimate the fraction of the "surfactant high" and "surfactant low" phenotype. The class label is given according to the phenotype with the highest estimate.

The function requires the vector `all_cpgs` of the CpG labels and a vector or matrix `data` of the corresponding beta values.

``` r
class_prediction <- surfactantClass(all_cpgs = sample$ID_REF, data = sample$VALUE)
```

    ## # NOTE - Selected GED algorithm: "lsfit" [signature data]
    ## Mapping signature ids onto target ids (method: auto) ... SKIP [signatures have no annotations]
    ## Limit/reorder to common set of features ... SKIP [missing feature names]
    ## Checking data dimension compatibility ... OK [4752 features x 2 cell types]
    ## Using cell type signatures:  [0 total]
    ## Checking log-scale ... data:YES - signatures:YES
    ## Reverting log-transform on signatures (base 2) ... OK
    ## Reverting log-transform on data (base 2) ... OK
    ## Normalizing signatures and target together (method: none) ... SKIP
    ##  Using ged algorithm: "lsfit"
    ##   Estimating cell proportions from cell-specific signatures [lsfit: ls]
    ##  Timing:
    ##    user  system elapsed 
    ##   1.753   0.018   1.799 
    ##  GED final wrap up ...  OK

``` r
class_prediction
```

    ##   low_estimate high_estimate           class
    ## 1    0.3457433     0.6542567 surfactant_high

The examplary LUAD sample can be classified to be "surfactant high".

#### Session info

``` r
sessionInfo()
```

    ## R version 3.5.2 (2018-12-20)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS Mojave 10.14.2
    ## 
    ## Matrix products: default
    ## BLAS: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.5/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_GB.UTF-8/en_GB.UTF-8/en_GB.UTF-8/C/en_GB.UTF-8/en_GB.UTF-8
    ## 
    ## attached base packages:
    ##  [1] stats4    compiler  parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] CellMix_1.6.2        GSEABase_1.42.0      graph_1.58.2        
    ##  [4] annotate_1.58.0      XML_3.98-1.16        AnnotationDbi_1.44.0
    ##  [7] IRanges_2.16.0       S4Vectors_0.20.1     stringr_1.4.0       
    ## [10] csSAM_1.2.4          NMF_0.22             Biobase_2.42.0      
    ## [13] BiocGenerics_0.28.0  cluster_2.0.7-1      rngtools_1.3.1      
    ## [16] pkgmaker_0.27        registry_0.5-1      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] matrixStats_0.54.0    bitops_1.0-6          bit64_0.9-7          
    ##  [4] doParallel_1.0.14     RColorBrewer_1.1-2    prabclus_2.2-7       
    ##  [7] tools_3.5.2           R6_2.4.0              DBI_1.0.0            
    ## [10] lazyeval_0.2.2        colorspace_1.4-1      trimcluster_0.1-2.1  
    ## [13] nnet_7.3-12           withr_2.1.2           tidyselect_0.2.5     
    ## [16] gridExtra_2.3         preprocessCore_1.42.0 bit_1.1-14           
    ## [19] diptest_0.75-7        scales_1.0.0          DEoptimR_1.0-8       
    ## [22] mvtnorm_1.0-10        robustbase_0.93-4     genefilter_1.62.0    
    ## [25] quadprog_1.5-5        digest_0.6.18         rmarkdown_1.11       
    ## [28] pkgconfig_2.0.2       htmltools_0.3.6       bibtex_0.4.2         
    ## [31] rlang_0.3.1           RSQLite_2.1.1         BiocInstaller_1.30.0 
    ## [34] mclust_5.4.3          gtools_3.8.1          dendextend_1.10.0    
    ## [37] dplyr_0.8.0.1         RCurl_1.95-4.11       magrittr_1.5         
    ## [40] modeltools_0.2-22     Matrix_1.2-15         Rcpp_1.0.1           
    ## [43] munsell_0.5.0         viridis_0.5.1         stringi_1.4.3        
    ## [46] whisker_0.3-2         yaml_2.2.0            MASS_7.3-51.1        
    ## [49] flexmix_2.3-15        plyr_1.8.4            grid_3.5.2           
    ## [52] blob_1.1.1            crayon_1.3.4          lattice_0.20-38      
    ## [55] splines_3.5.2         knitr_1.21            pillar_1.3.1         
    ## [58] fpc_2.1-11.1          corpcor_1.6.9         lpSolve_5.6.13       
    ## [61] reshape2_1.4.3        codetools_0.2-15      glue_1.3.1           
    ## [64] evaluate_0.12         foreach_1.4.4         gtable_0.3.0         
    ## [67] purrr_0.3.2           kernlab_0.9-27        assertthat_0.2.1     
    ## [70] ggplot2_3.1.0         xfun_0.4              gridBase_0.4-7       
    ## [73] limSolve_1.5.5.3      xtable_1.8-3          class_7.3-14         
    ## [76] survival_2.43-3       viridisLite_0.3.0     tibble_2.1.1         
    ## [79] iterators_1.0.10      beeswarm_0.2.3        memoise_1.1.0
