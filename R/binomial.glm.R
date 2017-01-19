
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

  # Take valid sample
  sample.fun     <- function(ydata, n, min.el = n/3) {
    ydata.factor <- factor(ydata)
    ydata.ix     <- seq_along(ydata)
    equal.list   <- list()
    equal.len    <- list()
    for(cl in levels(ydata.factor)) {
      equal.list[[cl]] <- ydata.ix[ydata.factor == cl]
      equal.len[[cl]]  <- length(equal.list[[cl]])
    }
    min.len.ix <- levels(ydata.factor)[sort(unlist(equal.len), index.return = T)$ix[1]]

    if (equal.len[[min.len.ix]] < min.el) {
      warning(sprintf('Cannot get a sample with at least %d elements', min.el))
    }
    #
    out.ix <- c()
    for(cl in levels(ydata.factor)) {
      out.ix <- c(out.ix, sample(equal.list[[cl]], min.el))
    }
    out.ix <- c(out.ix, sample(ydata.ix[-out.ix], n - length(out.ix)))
    #
    return(sort(out.ix))
  }

  # Compare Threshold
  threshold.cmp.fun <- function(error, threshold) {
    return(error <= threshold)
  }

  # Squared Error
  error.fun <- function(ydata, ydata.predicted) {
    return((ydata - ydata.predicted)^2)
  }

  # prediction function
  coef.fun <- function(object, ...) {
    coef(object = object)
  }

  # prediction function
  predict.fun <- function(object, newx, ...) {
    predict(object = object, newdata = data.frame(newx), type = 'response')
  }

  # Using RMSE
  model.error.fun <- function(ydata, ydata.predicted) {
    mean(sqrt(error.fun(ydata, ydata.predicted)))
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
    # Get Coefficients from model
    coef = coef.fun,
    # Compare error with threshold
    threshold.cmp = threshold.cmp.fun,
    # Model Error function name
    model.error.type = 'RMSE',
    # Error function name
    error.type = 'Squared error'))
}

