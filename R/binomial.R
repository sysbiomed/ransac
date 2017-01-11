
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
  fit.fun <- function(xdata, ydata, mc.cores = 1) {
    # need to suppress wanrings
    suppressWarnings(cv.glmnet(xdata, ydata, alpha = 0, family = 'binomial', lambda = c(0,1e-7, 1e-5, 1e-3, 1e-1)))
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
