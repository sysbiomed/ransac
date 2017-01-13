#
#' RANSAC
#'
#' @param xdata co-variates (rows: observations | coluns: co-variates)
#' @param ydata response
#' @param n minimum number of observations that will be fitted in the model
#' @param threshold error threshold to determine if observation is inliner
#' @param good.fit.perct discard models that have less than this percentage of inliners
#' @param k number of times model is run
#' @param family type of model that should be used
#' @param mc.cores number of processes to process the iteration of ransac
#'
#' @return
#' @export
#'
#' @examples
ransac <- function(xdata, ydata, n, threshold, good.fit.perct,
                   k = 100,
                   family = 'binomial',
                   mc.cores = 1, ...) {
  #
  # retrieve family by name
  if (is.character(family)) {
    family.fun <- switch(family,
                         binomial = family.binomial(),
                         binomial.glm = family.binomial.glm())
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }
  #
  flog.debug(paste0('Starting ransac with:\n',
                    '  k: %d\n',
                    '  n: %d\n',
                    '  threshold: %f\n',
                    '  good.fit.perct: %.2f'),
             k, n, threshold, good.fit.perct)
  #
  # Using parallel computing to speed up
  error.array <- mclapply(seq(k), function(ix) {
    flog.debug('%d / %d iteration', ix, k)
    #
    # draw a random sample from the dataset with size n
    xdata.min.ix <- family.fun$sample(ydata, n)
    xdata.min <- xdata[xdata.min.ix, ]
    ydata.min <- ydata[xdata.min.ix]
    #
    flog.debug('  ydata table: %d / %d (class 0 / class 1)', sum(ydata.min == 0), sum(ydata.min == 1))
    #
    # fit the model using the random sample
    inliers.model <- family.fun$fit.model(xdata.min, ydata.min, mc.cores = 1, ...)
    flog.debug('  finished the fitting %d', ix)
    #
    # test all the observations outside the random dataset
    #  if they are inliners of the model
    test.inliners.ix <- seq(nrow(xdata))[-xdata.min.ix]
    # TODO: drop the sapply and do this with matrices
    res <- sapply(test.inliners.ix, function(test.ix) {
      my.prediction <- family.fun$predict(object = inliers.model, newx = t(xdata[test.ix,]))
      family.fun$error(my.prediction, ydata[test.ix]) < threshold
    })
    also.inliners.ix <- test.inliners.ix[res]
    flog.debug('  also inliners: %d + %d = %d / %d (%.2f)',
               length(also.inliners.ix),
               nrow(xdata.min),
               length(also.inliners.ix) + nrow(xdata.min),
               nrow(xdata),
               good.fit.perct * nrow(xdata))
    #
    # Check if the model is good
    if (length(also.inliners.ix) + nrow(xdata.min) >= nrow(xdata) * good.fit.perct) {
      # this implies that we may have found a good model
      #  now test how good it is
      xdata.ix     <- c(also.inliners.ix, xdata.min.ix)

      #
      #
      # Calculating loss with original model fitted with n observations

      # only inliers + maybe.inliers
      my.residuals                      <- family.fun$predict(object = inliers.model, newx = xdata[xdata.ix,])
      error.inliers.model.inliers.ydata <- family.fun$model.error(ydata[xdata.ix], my.residuals)

      # with all xdata
      my.residuals                    <- family.fun$predict(object = inliers.model, newx = xdata)
      error.inliers.model.all.ydata   <- family.fun$model.error(ydata, my.residuals)

      #
      #
      # Calculating loss with refitted model with inliers + maybe.inliers

      # Re-estimating model
      all.inliers.model <- family.fun$fit.model(xdata[xdata.ix,], ydata[xdata.ix], mc.cores = 1, ...)

      # only inliers + maybe.inliers
      my.residuals                    <- family.fun$predict(object = all.inliers.model, newx = xdata[xdata.ix,])
      error.all.inliers.model.inliers <- family.fun$model.error(ydata[xdata.ix], my.residuals)

      # with all xdata
      my.residuals                      <- family.fun$predict(object = all.inliers.model, newx = xdata)
      error.all.inliers.model.all.ydata <- family.fun$model.error(ydata, my.residuals)

      #
      flog.debug(paste0('  %d / %d good fit!! There are %d total inliers Loss:\n',
                        '   initial model with inliers only: %f\n',
                        '       initial model with all data: %f\n',
                        '  refitted model with inliers only: %f\n',
                        '      refitted model with all data: %f'),
                 ix, k,
                 length(xdata.ix),
                 error.inliers.model.inliers.ydata,
                 error.inliers.model.all.ydata,
                 error.all.inliers.model.inliers,
                 error.all.inliers.model.all.ydata)
      #
      return(list(inliers.model     = inliers.model,
                  all.inliers.model = all.inliers.model,
                  #
                  error.inliers.model.inliers.ydata = error.inliers.model.inliers.ydata,
                  error.inliers.model.all.ydata     = error.inliers.model.all.ydata,
                  #
                  error.all.inliers.model.inliers   = error.all.inliers.model.inliers,
                  error.all.inliers.model.all.ydata = error.all.inliers.model.all.ydata,
                  #
                  min.ix     = xdata.min.ix,
                  inliers.ix = xdata.ix))
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
  best.error.inliers.inliers.ydata <- Inf
  best.model.inliers.inliers.ydata <- NA
  #
  best.model.inliers.all.ydata <- NA
  best.error.inliers.all.ydata <- Inf
  #
  best.model.all.inliers.all <- NA
  best.error.all.inliers.all <- Inf
  #
  best.model.all.inliers.consensus <- NA
  best.error.all.inliers.consensus <- Inf
  #
  best.model.all.inliers.model.inliers <- NA
  best.error.all.inliers.model.inliers <- Inf
  #
  lapply(seq_along(error.array), function(ix) {
    el <- error.array[[ix]]
    #
    if (ix %% 100 == 0) {
      flog.debug('Step %d / %d', ix, k)
    }
    if ((is.list(el) && length(el) == 0) || is.na(el) || is.null(el)) {
      return(0)
    }
    # model from initial set of inliers
    #  chosen by best RMSE with inliers dataset
    res <- lower.loss(el$error.inliers.model.inliers.ydata, el$inliers.model,
                      best.error.inliers.inliers.ydata, best.model.inliers.inliers.ydata)
    if (!is.null(res)) {
      best.error.inliers.inliers.ydata <<- res$best.error
      best.model.inliers.inliers.ydata <<- res$best.model
    }

    # model from initial set of inliers
    #  chosen by best RMSE with complete dataset
    res <- lower.loss(el$error.inliers.model.all.ydata, el$inliers.model,
                      best.error.inliers.all.ydata, best.model.inliers.all.ydata)
    if (!is.null(res)) {
      best.error.inliers.all.ydata <<- res$best.error
      best.model.inliers.all.ydata <<- res$best.model
    }

    #
    # models with inliers + maybe.inliers

    # model from all inliers (inliers + maybe.inliers)
    #  chosen by best RMSE with complete dataset
    res <- lower.loss(el$error.inliers.model.all.ydata, el$all.inliers.model,
                      best.error.all.inliers.all, best.model.all.inliers.all)
    if (!is.null(res)) {
      best.error.all.inliers.all <<- res$best.error
      best.model.all.inliers.all <<- res$best.model
    }

    # model from all inliers (inliers + maybe.inliers)
    #  chosen by best RMSE with inliers + maybe.inliers dataset
    res <- lower.loss(el$error.all.inliers.model.inliers, el$all.inliers.model,
                      best.error.all.inliers.model.inliers, best.model.all.inliers.model.inliers)
    if (!is.null(res)) {
      best.error.all.inliers.model.inliers <<- res$best.error
      best.model.all.inliers.model.inliers <<- res$best.model
    }

    # model from all inliers (inliers + maybe.inliers)
    #  chosen by most all inliers and resolves ties with RMSE
    res <- higher.consensus(family.fun,
                            el$error.all.inliers.model.all.ydata, el$all.inliers.model,
                            best.error.all.inliers.consensus, best.model.all.inliers.consensus)
    if (!is.null(res)) {
      best.error.all.inliers.consensus <<- res$best.error
      best.model.all.inliers.consensus <<- res$best.model
    }

    # trims down size of cache and memory usage
    error.array[[ix]]$inliers.model <<- NULL
    error.array[[ix]]$all.inliers.model <<- NULL
    #
    #
    return(1)
  })
  return(list(best.model.inliers.inliers.ydata = best.model.inliers.inliers.ydata,
              best.error.inliers.inliers.ydata = best.error.inliers.inliers.ydata,
              #
              best.model.inliers.all.ydata = best.model.inliers.all.ydata,
              best.error.inliers.all.ydata = best.error.inliers.all.ydata,
              #
              best.model.all.inliers.all.ydata = best.model.all.inliers.all,
              best.error.all.inliers.all.ydata = best.error.all.inliers.all,
              #
              best.model.all.inliers.model.inliers = best.model.all.inliers.model.inliers,
              best.error.all.inliers.model.inliers = best.error.all.inliers.model.inliers,
              #
              best.model.all.inliers.consensus = best.model.all.inliers.consensus,
              best.error.all.inliers.consensus = best.error.all.inliers.consensus,
              #
              errors = error.array ))
}


#' Choose model with higher consensus set
#'
#' @param family
#' @param el.error
#' @param el.model
#' @param best.model
#' @param best.error
#'
#' @return
#'
#' @examples
higher.consensus <- function(family.fun, el.error, el.model, best.error, best.model) {
  if (length(best.model) == 1 && is.na(best.model)) {
    best.nobs <- 0
  } else {
    best.nobs <- family.fun$nobs(best.model)
  }
  el.nobs   <- family.fun$nobs(el.model)
  if (!is.null(el.error) && !is.null(el.model) &&
      (best.nobs < el.nobs || (best.nobs == el.nobs && el.error < best.error))) {
    return(list(best.model = el.model,
                best.error = el.error))
  } else {
    return(NULL)
  }
}

#' Choose model with lower loss
#'
#' @param el.error
#' @param el.model
#' @param best.model
#' @param best.error
#'
#' @return
#'
#' @examples
lower.loss <- function(el.error, el.model, best.error, best.model) {
  if (!is.null(el.error) && !is.null(el.model) && el.error < best.error) {
    return(list(best.model = el.model,
                best.error = el.error))
  } else {
    return(NULL)
  }
}
