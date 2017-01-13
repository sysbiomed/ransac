
#' Support glm logistic regression
#'
#' Uses the binomial distribution
#'
#' @return
#' @export
#'
#' @examples
ransac.binomial.glm <- function() {
  #
  #
  # Functions

  nobs.fun <- function(model) {
    return(length(model$y))
  }

  # take valid sample
  sample.fun <- function(ydata, n, total.nruns = 1000, min.el = 8) {
    min.ix.count <- 0
    min.ix       <- NULL
    for (nrun in seq(total.nruns)) {
      ydata.size <- nrow(ydata)
      if (is.null(ydata.size)) {
        ydata.size <- length(ydata)
      }
      ix                <- sample(seq(ydata.size), n)
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
  error.fun <- function(ydata.predicted, ydata) {
    return((ydata.predicted - ydata)^2)
  }

  # prediction function
  coef.fun <- function(object) {
    coef(object = object)
  }

  # prediction function
  predict.fun <- function(object, newx) {
    predict(object = object, newdata = data.frame(newx), type = 'response')
  }

  # Using RMSE
  model.error.fun <- function(ydata.predicted, ydata) {
    mean(sqrt(error.fun(ydata.predicted, ydata)))
  }

  # fitting model
  fit.fun <- function(xdata, ydata, ...) {
    # generate balanced folds for the cross-validation

    new.df <- data.frame(xdata, logit_class = ydata)
    my.model <- glm(logit_class ~. , data = new.df, family = binomial(link = 'logit'),
                    control = glm.control(maxit = 1000))
    return(my.model)
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
    sample = sample.fun,
    # get observations used in model
    nobs = nobs.fun,
    #
    coef = coef.fun))
}

