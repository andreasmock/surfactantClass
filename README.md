Predicting the surfactant class from 450k arrays
================
Andreas Mock and Kolja Pocher

#### Introduction

The following tutorial describes the prediction of the"surfactant high" and "surfactant low" classes from methylome data (450k arrays). The predictor requires the `CellMix` and the `repmis` R package.

``` r
library(CellMix)
library(repmis)
```

The predictor can be loaded by:

``` r
source_data("https://github.com/andreasmock/surfactantClass/raw/master/surfactantClass.RData")
```

    ## Downloading data from: https://github.com/andreasmock/surfactantClass/raw/master/surfactantClass.RData

    ## SHA-1 hash of the downloaded data file is:
    ## d60f2c1b5f03a2c5a1003b156cc55bdcc9352973

    ## [1] "w"               "cpgs"            "surfactantClass"

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
    ##   1.798   0.012   1.828 
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
    ##  [1] repmis_0.5           CellMix_1.6.2        GSEABase_1.42.0     
    ##  [4] graph_1.58.2         annotate_1.58.0      XML_3.98-1.16       
    ##  [7] AnnotationDbi_1.44.0 IRanges_2.16.0       S4Vectors_0.20.1    
    ## [10] stringr_1.4.0        csSAM_1.2.4          NMF_0.22            
    ## [13] Biobase_2.42.0       BiocGenerics_0.28.0  cluster_2.0.7-1     
    ## [16] rngtools_1.3.1       pkgmaker_0.27        registry_0.5-1      
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] matrixStats_0.54.0    bitops_1.0-6          bit64_0.9-7          
    ##  [4] httr_1.4.0            doParallel_1.0.14     RColorBrewer_1.1-2   
    ##  [7] prabclus_2.2-7        R.cache_0.13.0        tools_3.5.2          
    ## [10] R6_2.4.0              DBI_1.0.0             lazyeval_0.2.2       
    ## [13] colorspace_1.4-1      trimcluster_0.1-2.1   nnet_7.3-12          
    ## [16] withr_2.1.2           tidyselect_0.2.5      gridExtra_2.3        
    ## [19] curl_3.3              preprocessCore_1.42.0 bit_1.1-14           
    ## [22] diptest_0.75-7        scales_1.0.0          DEoptimR_1.0-8       
    ## [25] mvtnorm_1.0-10        robustbase_0.93-4     genefilter_1.62.0    
    ## [28] quadprog_1.5-5        digest_0.6.18         R.utils_2.8.0        
    ## [31] rmarkdown_1.11        pkgconfig_2.0.2       htmltools_0.3.6      
    ## [34] bibtex_0.4.2          rlang_0.3.1           RSQLite_2.1.1        
    ## [37] BiocInstaller_1.30.0  mclust_5.4.3          gtools_3.8.1         
    ## [40] R.oo_1.22.0           dendextend_1.10.0     dplyr_0.8.0.1        
    ## [43] RCurl_1.95-4.11       magrittr_1.5          modeltools_0.2-22    
    ## [46] Matrix_1.2-15         Rcpp_1.0.1            munsell_0.5.0        
    ## [49] viridis_0.5.1         R.methodsS3_1.7.1     stringi_1.4.3        
    ## [52] whisker_0.3-2         yaml_2.2.0            MASS_7.3-51.1        
    ## [55] flexmix_2.3-15        plyr_1.8.4            grid_3.5.2           
    ## [58] blob_1.1.1            crayon_1.3.4          lattice_0.20-38      
    ## [61] splines_3.5.2         knitr_1.21            pillar_1.3.1         
    ## [64] fpc_2.1-11.1          corpcor_1.6.9         lpSolve_5.6.13       
    ## [67] reshape2_1.4.3        codetools_0.2-15      glue_1.3.1           
    ## [70] evaluate_0.12         data.table_1.12.0     foreach_1.4.4        
    ## [73] gtable_0.3.0          purrr_0.3.2           kernlab_0.9-27       
    ## [76] assertthat_0.2.1      ggplot2_3.1.0         xfun_0.4             
    ## [79] gridBase_0.4-7        limSolve_1.5.5.3      xtable_1.8-3         
    ## [82] class_7.3-14          survival_2.43-3       viridisLite_0.3.0    
    ## [85] tibble_2.1.1          iterators_1.0.10      beeswarm_0.2.3       
    ## [88] memoise_1.1.0
