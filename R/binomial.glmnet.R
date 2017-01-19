


#' Support glmnet logistic regression
#'
#' Uses the binomial distribution
#'
#' @return
#' @export
#'
#' @examples
ransac.binomial.glmnet <- function() {
  #
  #
  # Functions

  # Number of observation in model
  nobs.fun <- function(model) {
    return(model$nobs)
  }

  # Compare Threshold
  threshold.cmp.fun <- function(error, threshold) {
    return(error <= threshold)
  }

  # Squared Error
  error.fun <- function(ydata, ydata.predicted) {
    return((ydata - ydata.predicted)^2)
  }

  # Get Coefficients function
  coef.fun <- function(object, lambda, ...) {
    coef(object = object, s = lambda)
  }

  # Prediction function
  predict.fun <- function(object, newx, lambda, ...) {
    predict(object = object, newx = newx, s = lambda, type = 'response')
  }

  # Using RMSE
  model.error.fun <- function(ydata, ydata.predicted) {
    mean(sqrt(error.fun(ydata, ydata.predicted)))
  }

  # Fitting model
  fit.fun <- function(xdata, ydata, lambda, alpha = 0, penalty.factor = array(1, ncol(xdata)),
                      ...) {
    # need more than one lambda to guarantee convergence
    lambda.v <- sort(c(1000,100, 10, c(50, 10, 5, 4, 3, 2, 1.5, 1) * lambda), decreasing = T)

    # need to suppress wanrings
    return(suppressWarnings(glmnet(xdata, ydata,
                                   alpha = alpha,
                                   family = 'binomial',
                                   lambda = lambda.v,
                                   standardize = F,
                                   penalty.factor = penalty.factor)))
    #
  }

  parent.family <- ransac.binomial.glm()
  #
  return(list(
    # Squared error
    error = parent.family$error,
    # prediction function
    predict = predict.fun,
    # Using RMSE
    model.error = parent.family$model.error,
    # fitting model
    fit.model = fit.fun,
    # sample function
    sample = parent.family$sample,
    # get observations used in model
    nobs = nobs.fun,
    # Get coefficients from model
    coef = coef.fun,
    # threshold comparison
    threshold.cmp = threshold.cmp.fun,
    # Model Error function name
    model.error.type = 'RMSE',
    # Error function name
    error.type = 'Squared error'))
}

