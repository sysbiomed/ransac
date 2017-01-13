SOUND RANSAC
================

-   [Functions](#functions)
-   [Load BRCA data with 2 co-variates](#load-brca-data-with-2-co-variates)
    -   [Show SD of features](#show-sd-of-features)
-   [Running with k = 5000](#running-with-k-5000)
    -   [Description of parameters](#description-of-parameters)
    -   [RANSAC with the original data](#ransac-with-the-original-data)
    -   [RAMSAC with perturbation in first 10 and last 10](#ramsac-with-perturbation-in-first-10-and-last-10)

Functions
---------

Setup logging

``` r
#
# setup logger
flog.layout(layout.simple.parallel)
```

    ## NULL

``` r
flog.appender(appender.tee('logger.txt'))
```

    ## NULL

``` r
flog.threshold(DEBUG)
```

    ## NULL

Load BRCA data with 2 co-variates
---------------------------------

``` r
#
if (!exists('xdata') || !exists('ydata')) {
  genes <- c("ENSG00000072778","ENSG00000235505")
  load('BRCA.data.norm.RData')
  xdata <- BRCA.data.norm[,genes]
  ydata <- data.Y
  colnames(ydata) <- 'logit_class'
  data.name <- 'BRCA (2 co-variates)'

}

flog.info(paste0('Description of \'%s\' dataset:\n', '        cases: %d\n', '   class == 1: %d\n',
              '   class == 0: %d\n', '---------------------\n', '  co-variates: %d'),
          data.name, nrow(xdata), sum(ydata == 1), sum(ydata == 0),  ncol(xdata))
```

    ## INFO [2017-01-13 21:48:10 4813] Description of 'BRCA (2 co-variates)' dataset:
    ##         cases: 226
    ##    class == 1: 113
    ##    class == 0: 113
    ## ---------------------
    ##   co-variates: 2

### Show SD of features

``` r
if (ncol(xdata) > 100) {
  xdata.melt <- melt(xdata[, sample(seq(ncol(xdata)), 50)])
  flog.info('xdata has more than 50 co-variates (%d), showing a random subset of 50.', ncol(xdata))
} else {
  xdata.melt <- melt(xdata)
}
colnames(xdata.melt) <- c('id', 'feature', 'value')

xdata.melt$feature       <- factor(xdata.melt$feature)
xdata.melt$feature.short <- xdata.melt$feature
levels(xdata.melt$feature.short) <- gsub('ENSG0+','', levels(xdata.melt$feature))
#
g <- ggplot(xdata.melt, aes(feature.short, value)) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, face = 'italic', hjust = 0, vjust = .5)) + xlab('Feature') + 
  geom_boxplot(notch = TRUE, varwidth = TRUE, colour = my.colors(1), outlier.colour = my.colors(5), outlier.shape = 1)
# g + scale_y_continuous(trans = 'log1p') + ylab('Value (log1p scale)')
g + ylab('Value (no scale)')
```

![](ransac_files/figure-markdown_github/standard_deviation-1.svg)

Running with k = 5000
---------------------

``` r
# seed to make this reproducible
RNGkind("L'Ecuyer-CMRG")

# Parameters
if (exists('force.k'))    { n <- force.n }                 else { n <- 30 }
if (exists('force.k'))    { k <- force.k }                 else { k <- 100 }
if (exists('force.seed')) { my.seed <- force.seed }        else { my.seed <- 209 }
if (exists('force.thre')) { threshold <- force.thre }      else { threshold <- .35^2 }
if (exists('force.good')) { good.fit.perct <- force.good } else { good.fit.perct <- .8 }
if (exists('force.fami')) { family <- force.fami }         else { family <- 'binomial.glm' }
if (exists('force.core')) { mc.cores <- force.core }       else { mc.cores <- 16 }
```

### Description of parameters

-   `k`: number of iterations
-   `n`: cases to be included in the model (this should be the minimum number of observations required to fit the model)
-   `threshold`: error threshold to keep case as inlier (`error` &lt; `threshold`)
-   `good.fit.perct`: percentage of inlier cases that are required to consider a valid model

``` r
flog.info(paste0('Description of parameters:\n', 
                 '                k: %d\n', '                n: %d\n',
                 '        threshold: %g\n', '  good.fit.perct : %g'), 
          k, n , threshold, good.fit.perct)
```

    ## INFO [2017-01-13 21:48:11 4813] Description of parameters:
    ##                 k: 100
    ##                 n: 30
    ##         threshold: 0.1225
    ##   good.fit.perct : 0.8

### RANSAC with the original data

``` r
set.seed(my.seed)
#
cache.name <- generate.cache.name('results.ransac')
if (file.exists(cache.name)) {
  load(cache.name)
  flog.info('Loading cached ransac from %s', cache.name)
} else {
  results.ransac <- ransac(xdata, ydata, n = n, threshold = threshold, 
                           good.fit.perct = good.fit.perct, 
                           k = k, family = family, mc.cores = mc.cores,
                           nlambda = 100, alpha = 0)
  flog.info('Saving cache to disk (%s)', cache.name)
  save(results.ransac, file = cache.name)
}
```

    ## DEBUG [2017-01-13 21:48:11 4813] Starting ransac with:
    ##   k: 100
    ##   n: 30
    ##   threshold: 0.122500
    ##   good.fit.perct: 0.80
    ## DEBUG [2017-01-13 21:48:15 4813] Finished running all iterations, consolidating...
    ## DEBUG [2017-01-13 21:48:15 4813] Step 100 / 100
    ## INFO [2017-01-13 21:48:15 4813] Saving cache to disk (cache/results.ransac-226-2-100-30-0.122500-0.800000-209-BRCA_(2_co-variates).RData)

``` r
print.ransac.results(results.ransac, xdata, ydata, family = family, name = 'Original data', nlambda = 100)
```

    ## INFO [2017-01-13 21:48:15 4813] Results for k = 100 iterations

![](ransac_files/figure-markdown_github/test-1.svg)![](ransac_files/figure-markdown_github/test-2.svg)

    ## INFO [2017-01-13 21:48:16 4813] Information on RANSAC and Baseline model
    ## INFO [2017-01-13 21:48:16 4813]   2 Co-Var.  221 Obs 2.212389e-02 RMSE in RANSAC
    ## INFO [2017-01-13 21:48:16 4813] ----------------- Baseline ----------
    ## INFO [2017-01-13 21:48:16 4813]   2 Co-Var. 226 Obs 4.111275e-02 RMSE in Baseline
    ## INFO [2017-01-13 21:48:16 4813] 
    ## INFO [2017-01-13 21:48:16 4813] 
    ## INFO [2017-01-13 21:48:16 4813] False Positive/Negative
    ## INFO [2017-01-13 21:48:16 4813]   1 / 4 -- RANSAC
    ## INFO [2017-01-13 21:48:16 4813]   4 / 4 -- Baseline
    ## INFO [2017-01-13 21:48:16 4813] 
    ## INFO [2017-01-13 21:48:16 4813] 
    ## INFO [2017-01-13 21:48:16 4813] Misclassifications index
    ## INFO [2017-01-13 21:48:16 4813] 
    ## INFO [2017-01-13 21:48:16 4813]   RANSAC
    ## INFO [2017-01-13 21:48:16 4813]     False Positive 157
    ## INFO [2017-01-13 21:48:16 4813]     False Negative 56, 66, 75, 85
    ## INFO [2017-01-13 21:48:16 4813] 
    ## INFO [2017-01-13 21:48:16 4813]   Baseline
    ## INFO [2017-01-13 21:48:16 4813]     False Positive 134, 157, 206, 220
    ## INFO [2017-01-13 21:48:16 4813]     False Negative 22, 56, 60, 85

### RAMSAC with perturbation in first 10 and last 10

``` r
set.seed(my.seed)
#
ydata.pert <- ydata
ix.marta <- c(23,26,29,35,41,44,54,79,84,91,98,116,129,136,154,155,161,183,187,218)
ydata.pert[ix.marta] <- (ydata.pert[ix.marta] - 1) * -1
len <- length(ydata)
#
#ydata.pert[sample(which(ydata == 1),10)] <- 0
#ydata.pert[sample(which(ydata == 0),10)] <- 1
#
cache.name <- generate.cache.name('results.ransac.pert')
if (file.exists(cache.name)) {
  load(cache.name)
} else {
  results.ransac.pert <- ransac(xdata, ydata.pert, n = n, threshold = threshold, 
                                good.fit.perct = good.fit.perct, k = k,
                                family = family, mc.cores = mc.cores,
                                nlambda = 100, alpha = 0)
  flog.info('Saving cache to disk (%s)', cache.name)
  save(results.ransac.pert, file = cache.name)
}
```

    ## DEBUG [2017-01-13 21:48:16 4813] Starting ransac with:
    ##   k: 100
    ##   n: 30
    ##   threshold: 0.122500
    ##   good.fit.perct: 0.80
    ## DEBUG [2017-01-13 21:48:20 4813] Finished running all iterations, consolidating...
    ## DEBUG [2017-01-13 21:48:20 4813] Step 100 / 100
    ## INFO [2017-01-13 21:48:20 4813] Saving cache to disk (cache/results.ransac.pert-226-2-100-30-0.122500-0.800000-209-BRCA_(2_co-variates).RData)

``` r
print.ransac.results(results.ransac.pert, xdata, ydata.pert, 
                     family = family, 
                     name = 'With perturbation, error: ydata_perturbation - ydata_predicted')
```

    ## INFO [2017-01-13 21:48:20 4813] Results for k = 100 iterations

![](ransac_files/figure-markdown_github/test_perturbation-1.svg)![](ransac_files/figure-markdown_github/test_perturbation-2.svg)

    ## INFO [2017-01-13 21:48:21 4813] Information on RANSAC and Baseline model
    ## INFO [2017-01-13 21:48:21 4813]   2 Co-Var.  199 Obs 1.184836e-01 RMSE in RANSAC
    ## INFO [2017-01-13 21:48:21 4813] ----------------- Baseline ----------
    ## INFO [2017-01-13 21:48:21 4813]   2 Co-Var. 226 Obs 2.167997e-01 RMSE in Baseline
    ## INFO [2017-01-13 21:48:21 4813] 
    ## INFO [2017-01-13 21:48:21 4813] 
    ## INFO [2017-01-13 21:48:21 4813] False Positive/Negative
    ## INFO [2017-01-13 21:48:21 4813]   15 / 12 -- RANSAC
    ## INFO [2017-01-13 21:48:21 4813]   14 / 13 -- Baseline
    ## INFO [2017-01-13 21:48:21 4813] 
    ## INFO [2017-01-13 21:48:21 4813] 
    ## INFO [2017-01-13 21:48:21 4813] Misclassifications index
    ## INFO [2017-01-13 21:48:21 4813] 
    ## INFO [2017-01-13 21:48:21 4813]   RANSAC
    ## INFO [2017-01-13 21:48:21 4813]     False Positive 23, 26, 29, 35, 41, 44, 54, 79, 84, 91, 98, 134, 157, 206, 220
    ## INFO [2017-01-13 21:48:21 4813]     False Negative 22, 56, 60, 116, 129, 136, 154, 155, 161, 183, 187, 218
    ## INFO [2017-01-13 21:48:21 4813] 
    ## INFO [2017-01-13 21:48:21 4813]   Baseline
    ## INFO [2017-01-13 21:48:21 4813]     False Positive 23, 26, 29, 35, 41, 44, 54, 79, 84, 91, 98, 134, 157, 206
    ## INFO [2017-01-13 21:48:21 4813]     False Negative 56, 60, 75, 85, 116, 129, 136, 154, 155, 161, 183, 187, 218

``` r
#
```

#### How well does it recover the original data

From the dataset with perturbation, i.e.:

-   Trained with dirty dataset
-   Tested agains original

``` r
print.ransac.results(results.ransac.pert, xdata, ydata, ydata.original = ydata.pert,
                     family = family, 
                     name = 'With perturbation, error: ydata_original_class - ydata_predicted')
```

    ## INFO [2017-01-13 21:48:21 4813] Results for k = 100 iterations

![](ransac_files/figure-markdown_github/recover_it-1.svg)![](ransac_files/figure-markdown_github/recover_it-2.svg)

    ## INFO [2017-01-13 21:48:22 4813] Information on RANSAC and Baseline model
    ## INFO [2017-01-13 21:48:22 4813]   2 Co-Var.  199 Obs 2.998803e-02 RMSE in RANSAC
    ## INFO [2017-01-13 21:48:22 4813] ----------------- Baseline ----------
    ## INFO [2017-01-13 21:48:22 4813]   2 Co-Var. 226 Obs 1.523588e-01 RMSE in Baseline
    ## INFO [2017-01-13 21:48:22 4813] 
    ## INFO [2017-01-13 21:48:22 4813] 
    ## INFO [2017-01-13 21:48:22 4813] False Positive/Negative
    ## INFO [2017-01-13 21:48:22 4813]   4 / 3 -- RANSAC
    ## INFO [2017-01-13 21:48:22 4813]   3 / 4 -- Baseline
    ## INFO [2017-01-13 21:48:22 4813] 
    ## INFO [2017-01-13 21:48:22 4813] 
    ## INFO [2017-01-13 21:48:22 4813] Misclassifications index
    ## INFO [2017-01-13 21:48:22 4813] 
    ## INFO [2017-01-13 21:48:22 4813]   RANSAC
    ## INFO [2017-01-13 21:48:22 4813]     False Positive 134, 157, 206, 220
    ## INFO [2017-01-13 21:48:22 4813]     False Negative 22, 56, 60
    ## INFO [2017-01-13 21:48:22 4813] 
    ## INFO [2017-01-13 21:48:22 4813]   Baseline
    ## INFO [2017-01-13 21:48:22 4813]     False Positive 134, 157, 206
    ## INFO [2017-01-13 21:48:22 4813]     False Negative 56, 60, 75, 85
