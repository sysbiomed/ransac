


#' Support glmnet logistic regression
#'
#' Uses the binomial distribution
#'
#' @return
#' @export
#'
#' @examples
ransac.binomial.glmnet <- function(auc = F, residuals = 'pearson') {
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

  # Get Coefficients function
  coef.fun <- function(object, lambda, ...) {
    coef(object = object, s = lambda)
  }

  # Prediction function
  predict.fun <- function(object, newx, lambda, ...) {
    predict(object = object, newx = newx, s = lambda, type = 'response')
  }

  # Using AUC
  model.error.auc.fun <- function(ydata.predicted, ydata) {
    roc.dat <- AUC::roc(ydata.predicted, labels = factor(ydata))
    return(1 - AUC::auc(roc.dat))
  }

  # Using RMSE
  model.error.fun <- function(ydata.predicted, ydata) {
    mean(sqrt(error.fun(ydata, ydata.predicted)))
  }

  # Fitting model
  fit.fun <- function(xdata, ydata, lambda, alpha = 0, penalty.factor = array(1, ncol(xdata)),
                      intercept = TRUE,
                      ...) {
    # need more than one lambda to guarantee convergence
    lambda.v <- find.lambda(lambda)

    # need to suppress wanrings
    return(suppressWarnings(glmnet(xdata, ydata,
                                   alpha = alpha,
                                   family = 'binomial',
                                   lambda = lambda.v,
                                   standardize = F,
                                   intercept= FALSE,
                                   penalty.factor = penalty.factor)))
    #
  }

  parent.family <- ransac.binomial.glm(auc, residuals)
  #
  return(list(
    # Squared error
    error = parent.family$error,
    # prediction function
    predict = predict.fun,
    # Using RMSE
    model.error = if (auc) model.error.auc.fun else model.error.fun,
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
    # model name
    model.name = sprintf('GLMNET (%s)', parent.family$model.error.type),
    # Model Error function name
    model.error.type = parent.family$model.error.type,
    # Error function name
    error.type = parent.family$error.type))
}

