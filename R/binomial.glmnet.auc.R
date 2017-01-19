


#' Support glmnet logistic regression
#'
#' Uses the binomial distribution
#'
#' @return
#' @export
#'
#' @examples
ransac.binomial.glmnet.auc <- function() {
  #
  #
  # Functions

  # Using AUC
  model.error.fun <- function(ydata.predicted, ydata) {
    roc.dat <- AUC::roc(ydata.predicted, labels = factor(ydata))
    return(1 - AUC::auc(roc.dat))
  }

  #
  parent.family <- ransac.binomial.glmnet()
  #
  return(list(
    # Squared error
    error = parent.family$error,
    # prediction function
    predict = parent.family$predict,
    # Using RMSE
    model.error = model.error.fun,
    # fitting model
    fit.model = parent.family$fit.model,
    # sample function
    sample = parent.family$sample,
    # get observations used in model
    nobs = parent.family$nobs,
    # Get coefficients from model
    coef = parent.family$coef,
    # threshold comparison
    threshold.cmp = parent.family$threshold.cmp,
    # Model Error function name
    model.error.type = 'AUC',
    # Error function name
    error.type = 'Squared error'))
}

