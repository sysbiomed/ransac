print.ransac.results <- function(results.ransac, xdata, ydata, family = 'binomial', name = '') {
  #
  # retrieve family by name
  if (is.character(family)) {
    family.fun <- switch(family,
                         binomial = family.binomial())
  } else {
    family.fun <- family
  }
  if (is.null(family.fun)) {
    stop('family is not well defined, see documentation')
  }
  #
  if (results.ransac$best.error == Inf) {
    flog.info('Could not find a good fit')
  } else {
    models <- list()
    models$ransac <- results.ransac$best.model
    models$glmnet <- family.fun$fit.model(xdata, ydata)
    #
    my.fits <- list()
    my.fits$ransac <- family.fun$predict(models$ransac, xdata)
    my.fits$glmnet <- family.fun$predict(models$glmnet, xdata)
    #
    my.rmse <- list()
    my.rmse$ransac <- family.fun$model.error(my.fits$ransac, ydata)
    my.rmse$glmnet <- family.fun$model.error(my.fits$glmnet, ydata)
    #
    my.df <- data.frame()
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$ransac), type = rep('ransac', nrow(xdata))))
    my.df <- rbind(my.df, data.frame(ix = seq(nrow(xdata)), val = (ydata - my.fits$glmnet), type = rep('glmnet', nrow(xdata))))
    my.df$type <- factor(my.df$type)
    colnames(my.df) <- c('ix', 'value', 'type')
    #
    my.df.class <- data.frame()
    ix.sorted <- sort(ydata, index.return = T)$ix
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$ransac[ix.sorted]), type = rep('ransac', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (my.fits$glmnet[ix.sorted]), type = rep('glmnet', nrow(xdata))))
    my.df.class <- rbind(my.df.class, data.frame(ix = seq(nrow(xdata)), val = (ydata[ix.sorted]),          type = rep('ydata',  nrow(xdata))))
    my.df.class$type <- factor(my.df.class$type)
    colnames(my.df.class) <- c('ix', 'value', 'type')
    #
    # misclassifications
    miscl <- list()
    miscl$ransac <- sum(abs(my.fits$ransac - ydata) > .5)
    miscl$glmnet <- sum(abs(my.fits$glmnet - ydata) > .5)
    title.str <- paste0(name, ' : ransac misclass. = ', miscl$ransac, ' : glmnet misclass. = ', miscl$glmnet)
    #
    g <- ggplot(data = my.df.class) + theme_minimal() + ggtitle(title.str) + ylab('Classification') + xlab('Model')
    g <- g + geom_point(aes(ix, value, color = type))
    print(g)
    #
    g <- ggplot(data = my.df) + theme_minimal() + ggtitle(title.str) + ylab('Error') + xlab('Model')
    g <- g + geom_quasirandom(aes(type, value, color = type))
    print(g)
    #
    flog.info('RMSE:\n   ransac = %f\n   glmnet = %f', my.rmse$ransac, my.rmse$glmnet)
  }
}
