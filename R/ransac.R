#
#' RANSAC
#'
#' @param xdata co-variates (rows: observations | coluns: co-variates)
#' @param ydata response
#' @param n minimum number of observations that will be fitted in the model
#' @param threshold error threshold to determine if observation is inliner
#' @param good.fit.perct discard models that have less than this percentage of inliners
#' @param k number of times model is run
#' @param family
#' @param mc.cores
#'
#' @return
#' @export
#'
#' @examples
ransac <- function(xdata, ydata, n, threshold, good.fit.perct,
                   k = 100,
                   family = 'binomial',
                   mc.cores = 1) {
  #
  # retrieve family by name
  if (is.character(family)) {
    family.fun <- switch(family,
                         binomial = family.binomial())
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }
  #
  flog.debug('Starting ransac with n = %d, threshold = %f and good.fit.perct = %f...',
             n, threshold, good.fit.perct)
  #
  # Using parallel computing to speed up
  error.array <- mclapply(seq(k), function(ix) {
    flog.debug('Running iteration %d / %d', ix, k)
    #
    # draw a random sample from the dataset with size n
    xdata.min.ix <- family.fun$sample(ydata, n)
    xdata.min <- xdata[xdata.min.ix, ]
    ydata.min <- ydata[xdata.min.ix]
    #
    flog.debug('  ydata table: %d / %d (class 0 / class 1)', sum(ydata.min == 0), sum(ydata.min == 1))
    #
    # fit the model using the random sample
    my.fit <- family.fun$fit.model(xdata.min, ydata.min)
    flog.debug('  finished the fitting %d', ix)
    #
    # test all the observations outside the random dataset
    #  if they are inliners of the model
    test.inliners.ix <- seq(nrow(xdata))[-xdata.min.ix]
    # TODO: drop the sapply and do this with matrices
    res <- sapply(test.inliners.ix, function(test.ix) {
      my.prediction <- family.fun$predict(object = my.fit, newx = t(xdata[test.ix,]))
      family.fun$error(my.prediction, ydata[test.ix]) < threshold
    })
    also.inliners.ix <- test.inliners.ix[res]
    flog.debug('  also inliners: %d + %d = %d / %d',
               length(also.inliners.ix),
               nrow(xdata.min),
               length(also.inliners.ix) + nrow(xdata.min),
               nrow(xdata))
    #
    # Check if the model is good
    if (length(also.inliners.ix) + nrow(xdata.min) >= nrow(xdata) * good.fit.perct) {
      # this implies that we may have found a good model
      #  now test how good it is
      xdata.ix     <- c(also.inliners.ix, xdata.min.ix)
      my.residuals <- family.fun$predict(object = my.fit, newx = xdata[xdata.ix,])
      this.error   <- family.fun$model.error(ydata[xdata.ix], my.residuals)
      flog.debug('  good fit!! error = %f', this.error)
      return(list(this.error = this.error, this.model = my.fit))
    }
    flog.debug('  did not have a good fit')
    return(list())
  }
  # parameters for mclapply
  , mc.cores = mc.cores, mc.allow.recursive = F, mc.cleanup = T,
  mc.set.seed = T)
  flog.debug('Finished running all iterations, consolidating...')
  #
  # Consolidate results finding only the best one
  #
  # Setting temporary variables to keep the best model
  best.error <- Inf
  best.model <- NA
  #
  lapply(seq_along(error.array), function(ix) {
    el <- error.array[[ix]]
    if (!is.null(el) && !is.na(el) && length(el) > 1 &&
        !is.null(el$this.error) && el$this.error < best.error) {
      best.error <<- el$this.error
      best.model <<- el$this.model
    }
  })
  return(list(best.model = best.model, best.error = best.error, errors = error.array ))
}
