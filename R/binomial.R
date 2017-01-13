
family.binomial <- function() {
  #
  #
  # Functions

  # take valid sample
  sample.fun = function(ydata, n, total.nruns = 1000, min.el = 8) {
    min.ix.count <- 0
    min.ix       <- NULL
    for (nrun in seq(total.nruns)) {
      ix                <- sample(seq(length(ydata)), n)
      count.table       <- table(ydata[ix])
      min.ix.count.temp <- min(count.table)
      if (length(count.table) >= 2){
        if (min.ix.count.temp >= min.el) {
          return(ix)
        } else if (min.ix.count.temp > min.ix.count) {
          min.ix       <- ix
          min.ix.count <- min.ix.count.temp
        }
      }
    }
    if (min(table(ydata[ix])) < 5) {
      stop('Could not find a good sample from dataset after.')
    }
    flog.warn('One class does not have minimum of 8 samples, it has %d', min.ix.count)
    return(min.ix)
  }

  # Squared Error
  error.fun = function(ydata.predicted, ydata) {
    return((ydata.predicted - ydata)^2)
  }

  # prediction function
  predict.fun = function(object, newx) {
    predict(object = object, newx = newx, s = 'lambda.min', type = 'response')
  }

  # Using RMSE
  model.error.fun = function(ydata.predicted, ydata) {
    mean(sqrt(error.fun(ydata.predicted, ydata)))
  }

  # fitting model
  fit.fun <- function(xdata, ydata, alpha = 0, nlambda = NULL, mc.cores = 1) {
    # generate balanced folds for the cross-validation

    # get index of each class
    ydata.class.0 <- which(ydata == levels(factor(ydata))[1])
    ydata.class.1 <- which(ydata == levels(factor(ydata))[2])
    # calculate folds for each class
    out.balanced <- balanced.cv.folds(ydata.class.0, ydata.class.1, 10)
    # join them
    foldid <- array(0,length(ydata))
    foldid[ydata.class.0] <- out.balanced$output[[1]]
    foldid[ydata.class.1] <- out.balanced$output[[2]]
    # need to suppress wanrings
    if (is.null(nlambda)) {
      # little regularization
      if (any(names(formals(get("cv.glmnet"))) == 'mc.cores')){
        # https://github.com/averissimo/rpackage-glmnet
        #  that uses parallel package instead
        return(suppressWarnings(cv.glmnet(xdata, ydata, alpha = alpha,
                                          foldid = foldid,
                                          family = 'binomial',
                                          lambda = c(1e-9, 1e-7, 1e-5, 1e-3, 1e-1),
                                          mc.cores = mc.cores)))
      } else {
        # CRAN's glmnet package
        return(suppressWarnings(cv.glmnet(xdata, ydata, alpha = alpha,
                                          foldid = foldid,
                                          family = 'binomial',
                                          lambda = c(0,1e-7, 1e-5, 1e-3, 1e-1))))
      }
    } else {
      # normal regularization
      if (any(names(formals(get("cv.glmnet"))) == 'mc.cores')){
        return(suppressWarnings(cv.glmnet(xdata, ydata, alpha = alpha,
                                          foldid = foldid,
                                          family = 'binomial',
                                          nlambda = nlambda,
                                          mc.cores = mc.cores)))
      } else {
        # CRAN's glmnet package
        return(suppressWarnings(cv.glmnet(xdata, ydata, alpha = alpha,
                                          foldid = foldid,
                                          family = 'binomial',
                                          nlambda = nlambda)))
      }
    }
    #
  }

  #
  #
  return(list(
    # Squared error
    error = error.fun,
    # prediction function
    predict = predict.fun,
    # Using RMSE
    model.error = model.error.fun,
    # fitting model
    fit.model = fit.fun,
    # sample function
    sample = sample.fun))
}

