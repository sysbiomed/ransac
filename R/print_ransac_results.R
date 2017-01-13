#' Function to print the results from RANSAC
#'
#' @param results.ransac
#' @param xdata
#' @param ydata
#' @param family
#' @param name
#' @param ydata.original
#' @param show.title
#' @param only_consensus
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot.ransac <- function(results.ransac, xdata, ydata,
                        family = 'binomial', name = '', ydata.original = NULL, show.title = T,
                        only_consensus = T, ...) {
  #
  # retrieve family by name
  if (is.character(family)) {
    family.fun <- switch(family,
                         binomial = ransac.binomial(),
                         binomial.glm = ransac.binomial.glm())
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }
  #
  flog.info('Results for k = %d iterations', length(results.ransac$errors))
  #
  models <- list()
  if (!only_consensus) {
    models$ransac <- results.ransac$best.model.inliers.inliers.ydata
    models$rans.2 <- results.ransac$best.model.inliers.all.ydata
    models$rans.3 <- results.ransac$best.model.all.inliers.model.inliers
    models$rans.4 <- results.ransac$best.model.all.inliers.all.ydata
  }
  models$rans.5 <- results.ransac$best.model.all.inliers.consensus
  if(is.null(ydata.original)) {
    models$glmnet <- family.fun$fit.model(xdata, ydata, ...)
  } else {
    models$glmnet <- family.fun$fit.model(xdata, ydata.original, ...)
  }
  #
  my.fits <- list()
  if (!only_consensus) {
    my.fits$ransac <- family.fun$predict(models$ransac, xdata)
    my.fits$rans.2 <- family.fun$predict(models$rans.2, xdata)
    my.fits$rans.3 <- family.fun$predict(models$rans.3, xdata)
    my.fits$rans.4 <- family.fun$predict(models$rans.4, xdata)
  }
  my.fits$rans.5 <- family.fun$predict(models$rans.5, xdata)
  my.fits$glmnet <- family.fun$predict(models$glmnet, xdata)
  #
  my.rmse <- list()
  if (!only_consensus) {
    my.rmse$ransac <- family.fun$model.error(my.fits$ransac, ydata)
    my.rmse$rans.2 <- family.fun$model.error(my.fits$rans.2, ydata)
    my.rmse$rans.3 <- family.fun$model.error(my.fits$rans.3, ydata)
    my.rmse$rans.4 <- family.fun$model.error(my.fits$rans.4, ydata)
  }
  my.rmse$rans.5 <- family.fun$model.error(my.fits$rans.5, ydata)
  my.rmse$glmnet <- family.fun$model.error(my.fits$glmnet, ydata)
  #
  description <- list()
  if (!only_consensus) {
    description$ransac <- 'RANSAC Initial + Inliers'
    description$rans.2 <- 'RANSAC Initial + {In,Out}liers'
    description$rans.3 <- 'RANSAC Refitted + Inliers'
    description$rans.4 <- 'RANSAC Refitted + {In,Out}liers'
    description$rans.5 <- 'RANSAC Refitted + Consensus'
  } else {
    description$rans.5 <- 'RANSAC'
  }
  description$glmnet <- 'Baseline'
  #
  my.df <- data.frame()
  if (!only_consensus) {
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$ransac), type = rep(description$ransac, nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$rans.2), type = rep(description$rans.2, nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$rans.3), type = rep(description$rans.3, nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$rans.4), type = rep(description$rans.4, nrow(xdata))))
  }
  my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$rans.5), type = rep(description$rans.5, nrow(xdata))))
  my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$glmnet), type = rep(description$glmnet, nrow(xdata))))
  my.df$type <- factor(my.df$type)
  colnames(my.df) <- c('ix', 'value', 'type')
  #
  my.df.class <- data.frame()
  ix.sorted <- sort(ydata, index.return = T)$ix
  if (!only_consensus) {
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$ransac[ix.sorted]), type = rep(description$ransac, nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$rans.2[ix.sorted]), type = rep(description$rans.2, nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$rans.3[ix.sorted]), type = rep(description$rans.3, nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$rans.4[ix.sorted]), type = rep(description$rans.4, nrow(xdata))))
  }
  my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$rans.5[ix.sorted]), type = rep(description$rans.5, nrow(xdata))))
  my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$glmnet[ix.sorted]), type = rep(description$glmnet, nrow(xdata))))
  my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (ydata[ix.sorted]),          type = rep('Observed',  nrow(xdata))))
  my.df.class$type <- factor(my.df.class$type)
  colnames(my.df.class) <- c('ix', 'value', 'type')
  #
  # misclassifications
  miscl <- list()
  if (!only_consensus) {
    miscl$ransac <- list()
    miscl$rans.2 <- list()
    miscl$rans.3 <- list()
    miscl$rans.4 <- list()
  }
  miscl$rans.5 <- list()
  miscl$glmnet <- list()
  #
  if (!only_consensus) {
    miscl$ransac$false.pos <- round(ydata - my.fits$ransac) == -1
    miscl$ransac$false.neg <- round(ydata - my.fits$ransac) == 1
    #
    miscl$rans.2$false.pos <- round(ydata - my.fits$rans.2) == -1
    miscl$rans.2$false.neg <- round(ydata - my.fits$rans.2) == 1
    #
    miscl$rans.3$false.pos <- round(ydata - my.fits$rans.3) == -1
    miscl$rans.3$false.neg <- round(ydata - my.fits$rans.3) == 1
    #
    miscl$rans.4$false.pos <- round(ydata - my.fits$rans.4) == -1
    miscl$rans.4$false.neg <- round(ydata - my.fits$rans.4) == 1
  }
  #
  miscl$rans.5$false.pos <- round(ydata - my.fits$rans.5) == -1
  miscl$rans.5$false.neg <- round(ydata - my.fits$rans.5) == 1
  #
  miscl$glmnet$false.pos <- round(ydata - my.fits$glmnet) == -1
  miscl$glmnet$false.neg <- round(ydata - my.fits$glmnet) == 1
  #
  g <- ggplot(data = my.df.class) + theme_minimal() + ylab('Classification') + xlab('Model')
  if (show.title) { g <- g + ggtitle(name) }
  g <- g + geom_point(aes(ix, value, color = type))
  print(g)
  #
  g <- ggplot(data = my.df) + theme_minimal() + ylab('Error') + xlab('Model')
  if (show.title) { g <- g + ggtitle(name) }
  g <- g + geom_quasirandom(aes(type, value, color = type), bandwidth = 2, method = 'quasirandom') + scale_y_continuous(limits = c(-1,1))
  print(g)
  #
  if (!only_consensus) {
    non.zero.ransac <- family.fun$coef(models$ransac)[-1] != 0
    non.zero.rans.2 <- family.fun$coef(models$rans.2)[-1] != 0
    non.zero.rans.3 <- family.fun$coef(models$rans.3)[-1] != 0
    non.zero.rans.4 <- family.fun$coef(models$rans.4)[-1] != 0
  }
  non.zero.rans.5 <- family.fun$coef(models$rans.5)[-1] != 0
  non.zero.glmnet <- family.fun$coef(models$glmnet)[-1] != 0
  #
  flog.info('Information on RANSAC and Baseline model')
  if (!only_consensus) {
    flog.info('  %d Co-Var. % 4d Obs %e RMSE in %s', sum(non.zero.ransac), family.fun$nobs(models$ransac), my.rmse$ransac, description$ransac)
    flog.info('  %d Co-Var. % 4d Obs %e RMSE in %s', sum(non.zero.rans.2), family.fun$nobs(models$rans.2), my.rmse$rans.2, description$rans.2)
    flog.info('  %d Co-Var. % 4d Obs %e RMSE in %s', sum(non.zero.rans.2), family.fun$nobs(models$rans.2), my.rmse$rans.2, description$rans.2)
    flog.info('  %d Co-Var. % 4d Obs %e RMSE in %s', sum(non.zero.rans.4), family.fun$nobs(models$rans.4), my.rmse$rans.4, description$rans.4)
  }
  flog.info('  %d Co-Var. % 4d Obs %e RMSE in %s', sum(non.zero.rans.5), family.fun$nobs(models$rans.5), my.rmse$rans.5, description$rans.5)
  flog.info('----------------- Baseline ----------')
  flog.info('  %d Co-Var. %d Obs %e RMSE in Baseline', sum(non.zero.glmnet), family.fun$nobs(models$glmnet), my.rmse$glmnet)
  #
  flog.info('')
  flog.info('')
  flog.info('False Positive/Negative')
  if (!only_consensus) {
    flog.info('  %d / %d -- %s', sum(miscl$ransac$false.pos), sum(miscl$ransac$false.neg), description$ransac)
    flog.info('  %d / %d -- %s', sum(miscl$rans.2$false.pos), sum(miscl$rans.2$false.neg), description$rans.2)
    flog.info('  %d / %d -- %s', sum(miscl$rans.3$false.pos), sum(miscl$rans.3$false.neg), description$rans.3)
    flog.info('  %d / %d -- %s', sum(miscl$rans.4$false.pos), sum(miscl$rans.4$false.neg), description$rans.4)
  }
  flog.info('  %d / %d -- %s', sum(miscl$rans.5$false.pos), sum(miscl$rans.5$false.neg), description$rans.5)
  flog.info('  %d / %d -- %s', sum(miscl$glmnet$false.pos), sum(miscl$glmnet$false.neg), description$glmnet)
  #
  flog.info('')
  flog.info('')
  flog.info('Misclassifications index')
  if (!only_consensus) {
    flog.info('')
    flog.info('  %s', description$ransac)
    flog.info('    False Positive %s', paste(which(miscl$ransac$false.pos), collapse = ', '))
    flog.info('    False Negative %s', paste(which(miscl$ransac$false.neg), collapse = ', '))
    #
    flog.info('')
    flog.info('  %s', description$rans.2)
    flog.info('    False Positive %s', paste(which(miscl$rans.2$false.pos), collapse = ', '))
    flog.info('    False Negative %s', paste(which(miscl$rans.2$false.neg), collapse = ', '))
    #
    flog.info('')
    flog.info('  %s', description$rans.3)
    flog.info('    False Positive %s', paste(which(miscl$rans.3$false.pos), collapse = ', '))
    flog.info('    False Negative %s', paste(which(miscl$rans.3$false.neg), collapse = ', '))
    #
    flog.info('')
    flog.info('  %s', description$rans.4)
    flog.info('    False Positive %s', paste(which(miscl$rans.4$false.pos), collapse = ', '))
    flog.info('    False Negative %s', paste(which(miscl$rans.4$false.neg), collapse = ', '))
  }
  #
  flog.info('')
  flog.info('  %s', description$rans.5)
  flog.info('    False Positive %s', paste(which(miscl$rans.5$false.pos), collapse = ', '))
  flog.info('    False Negative %s', paste(which(miscl$rans.5$false.neg), collapse = ', '))
  #
  flog.info('')
  flog.info('  %s', description$glmnet)
  flog.info('    False Positive %s', paste(which(miscl$glmnet$false.pos), collapse = ', '))
  flog.info('    False Negative %s', paste(which(miscl$glmnet$false.neg), collapse = ', '))

}

