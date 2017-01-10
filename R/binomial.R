
family.binomial <- function() {
  #
  #
  # Functions

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
    fit.model = fit.fun))
}
