
#' Support glm logistic regression
#'
#' Uses the binomial distribution
#'
#' @return
#' @export
#'
#' @examples
ransac.binomial.glm <- function(auc = F, residuals = 'pearson') {
  #
  #
  # Functions

  nobs.fun <- function(model) {
    return(length(model$y))
  }

  # Take valid sample
  sample.fun     <- function(ydata, n, min.el = max(c(1,round(n/3)))) {
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

  # Pearson
  error.dev.fun <- function(ydata, ydata.predicted) {
    top <- (ydata - ydata.predicted)
    top.ix.0 <- top == 0
    predicted.0.or.1 <- ydata.predicted == 1 | ydata.predicted == 0
    result <- abs(top) / sqrt(ydata.predicted*(1-ydata.predicted))
    result[top.ix.0] <- 0
    result[predicted.0.or.1] <- Inf
    return(result)
  }

  # prediction function
  coef.fun <- function(object, ...) {
    coef(object = object)
  }

  # prediction function
  predict.fun <- function(object, newx, ...) {
    predict(object = object, newdata = data.frame(newx), type = 'response')
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

  # fitting model
  fit.fun <- function(xdata, ydata, intercept = T, ...) {
    # generate balanced folds for the cross-validation

    new.df <- data.frame(xdata, logit_class = ydata)
    suppressWarnings({
      my.model <- glm(logit_class ~. , data = new.df, family = binomial(link = 'logit'),
                      control = glm.control(maxit = 1000),
                      intercept = intercept)
    })
    my.model$data <- NULL
    return(my.model)
    #
  }

  model.error.type <- if (auc) 'AUC' else 'Squared Error'
  #
  error.type.fun <- switch(residuals,
                           squared.error = error.fun,
                           pearson       = error.dev.fun)
  error.type.txt <- switch(residuals,
                           squared.error = 'Squared Error',
                           pearson       = 'Pearson Residuals')

  #
  return(list(
    # Squared error
    error = error.type.fun,
    # prediction function
    predict = predict.fun,
    # Using RMSE
    model.error = if (auc) model.error.auc.fun else model.error.fun,
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
    # model name
    model.name = sprintf('GLM (%s)', model.error.type),
    # Model Error function name
    model.error.type = model.error.type,
    # Error function name
    error.type = error.type.txt))
}

